#ifndef MLFMM_Env_FMM_CPP
#define MLFMM_Env_FMM_CPP

#include "MLFMM_Env.h"

namespace FMM {


  Indexing MLFMM_Env::global_index(LMAX);
  

  MLFMM_Env::MLFMM_Env( HelmKernel K_, complex k_out_, std::vector<Vec3> p,
			std::vector<Scatterer> ScL, int levels_, double eps)
    : K(K_), k_out(k_out_), leaf_Translation(NScat), sourcefield(p), 
      tree( p, levels_ ),    // Create the tree
      levels( levels_+1 ),       // Initialize std::vector of levels
      Interface_Level(-1)
  {
    
    if( verbose >= 2 )
      cerr << tree << endl;
    
    
    int nLevels = tree.numLevels();                   // Number of levels in the tree
    //double cutoff =  CUTOFF_SIZE*2*PI/real(k_out);    // Interface between low- and high-frequency
    //cout << "cutoff : " << cutoff << endl;
    
    // Construct each Level and give it the boxes it owns
    for( int L = 2; L <= nLevels; ++L ) {
      levels[L] = Level(tree.getBoxSize(L), k_out, CUTOFF_SIZE);
      for( int n = 0; n < (1 << (DIM*L)); ++n ) {
	Box* b = getBox(n,L);
	if( b != NULL ) getLevel(L).addBox(b);
      }
    }  
    if( verbose >= 2 ) cerr << "   Building Transfer and Close Lists..." << endl;
    
    for( int L = 2; L <= nLevels; ++L )
      getLevel(L).setL();
    
    // Recursively build transfer and close lists from rootbox
    defineTransferAndClose(getBox(0,0), 0); 
    if( verbose >= 2 ) cerr << "   Transfer and Close Lists Built" << endl << endl;
    
    
    //TODO : Should only be done for levels above cutoff line
    // Compute the quadrature for each level of the tree and Initialize fields
    getLevel(nLevels).LEAF = true;

    
    for( int L = 2; L <= nLevels; ++L ) {
      getLevel(L).defineQuadrature( K, eps );
      
      // Establish which level is at the interface between low- and high-frequency
	    if( Interface(getLevel(L).getBoxSize()/2, CUTOFF_SIZE) || (Interface_Level == -1 && L == nLevels && getLevel(L-1).HF ) )
	      Interface_Level = L;
	    
	    getLevel(L).initializeFields();
    }
    
    if( Interface_Level != -1 ) getLevel(Interface_Level).INTF = true;     // Keep track of interface level;
    
    if( verbose >= 2 ){
      cout << "      Interface level : " << Interface_Level << endl;
      cout << "   Quadratures defined" << endl;
    }
    
    
    // Compute SHT at interface level
    if( Interface_Level != nLevels && Interface_Level >= 0 )
      getLevel( Interface_Level ).defineSH_T( getLevel( Interface_Level+1 ).getQuad(), getLevel( Interface_Level+1 ).getIndex());
    else if( Interface_Level >= 0 )
      leaf_SH_T = new SH_Transform(getLevel(nLevels).getQuad(), &global_index);
    
    if( verbose >= 2 ) cout << "   SHT defined" << endl;
    
    
    // TODO: NEED FIX HERE
    // Initialize S_Rotation in Transfer Utilities
    int maxIndexSize = LMAX;
    for( int L = 2; L <= nLevels; ++L ){
      if( !getLevel(L).HF || L == Interface_Level ){
	cout << (*getLevel(L).getIndex())(getLevel(L).getIndex()->size()-1,0) << endl;
	maxIndexSize = std::max(maxIndexSize, (*getLevel(L).getIndex())(getLevel(L).getIndex()->size()-1,0) );
      }
    }
    
    
    // TODO: Currently making it bigger than it needs to be?
    // Initialiaze static variables for Transfer Utils
    S_Rotation::S_Init( 2*maxIndexSize );
    Z_transfer::Z_Init( 2*maxIndexSize );
    if( verbose >= 2 ) cerr << "   S_Rotation defined" << endl;
    
    // Define the interpolators and anterpolators required for each level
    for( int L = 2; L <= nLevels; ++L ) {	  
      if( L != 2 ) 
	getLevel(L).defineInterp( getLevel(L-1).getQuad() );
      if( L != nLevels )
	getLevel(L).defineAnterp( getLevel(L+1).getQuad() );
    }
    
    if( verbose >= 2 )
      cout << "   Reterpolators Defined" << endl;
    
    
    
    // Compute the transfer functions for each level of the tree
    for( int L = 2; L <= nLevels; ++L ) {
      cout << L << endl;
      getLevel(L).defineTransfers( eps, k_out );
    }
    
    if( verbose >= 2 )
      cout << "   Transfer Functions Defined" << endl;
    
    
    // Compute the translation operators for each level of the tree 
    Indexing *pIndex, *cIndex;
    for( int L = 2; L < nLevels; ++L ) {
      (L==2) ? pIndex = NULL : pIndex = getLevel(L-1).getIndex();
      cIndex = getLevel(L+1).getIndex();
      
      getLevel(L).defineTranslations(pIndex, cIndex, k_out);
    }
    
    if( verbose >= 2 )
      cout << "   Translation Functions (at higher levels) Defined" << endl;
    
    
    
    // Compute the translation operators for each scatterer at the leaf level
    // TODO : Should probably be in class 'Level'
    // For each box at leaf level

    //Indexing *leafIndex; leafIndex = new Indexing(LMAX_VSWF);
    Level& leafLevel = getLevel(nLevels);
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;
      const Vec3& center = b->center;
      const Box::pointIter piEnd = b->pointIndex.end();
      
      // For each particle in box
      for( Box::pointIter pi = b->pointIndex.begin(); pi != piEnd; ++pi ) {
	int idx = *pi;
	const Vec3 r = center - sourcefield[idx];
	
	if( Interface_Level == nLevels )
	  leaf_Translation[idx] = new Leaf_HF_Translation_Function(r, k_out, getLevel(nLevels).getQuad());
	else
	  leaf_Translation[idx] = new Leaf_Direct_Translation_Function(r, k_out, getLevel(nLevels).getIndex());

	//cout << "leaf level size: " <<  (getLevel(nLevels).getIndex())->size() << endl;
	
      }
    }
    
    
    
    if( verbose >= 2 )
      cout << "   Translation Functions (leaf level) Defined" << endl;
    
    
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;
      
      // for all interior points
      const Box::pointIter nEnd = b->pointIndex.end();
      for( Box::pointIter ni = b->pointIndex.begin(); ni != nEnd; ++ni ) {
	int n = *ni;
	
	// Add close transfers from other points inside this box
	for( Box::pointIter mi = b->pointIndex.begin(); mi != ni; ++mi ) {
	  int m = *mi;
	  
	  //from m to n
	  //double begin = clock();
	  Vec3 r0 = sourcefield[n] - sourcefield[m];
	  leafLevel.addCloseTransfer(n, m, r0, k_out);
	  //cout << "close transfer construct : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
	  
	  //from n to m
	  Vec3 r1 = -r0;
	  leafLevel.addCloseTransfer(m, n, r1, k_out);
	}
      }
    }
    
    
    // Add close transfers from other points inside neighboring boxes
    for( CVIter cvi = leafLevel.closeList.begin(); cvi != leafLevel.closeList.end(); ++cvi ) {
      Box* b1 = cvi->b1;
      Box* b2 = cvi->b2;
      
      // For all pairs of points inside box b1 and b2
      const Box::pointIter nStart = b1->pointIndex.begin();
      const Box::pointIter nEnd   = b1->pointIndex.end();
      const Box::pointIter mStart = b2->pointIndex.begin();
      const Box::pointIter mEnd   = b2->pointIndex.end();
      
      for( Box::pointIter ni = nStart; ni != nEnd; ++ni ) {
	int n = *ni;
	const Vec3& pn = sourcefield[n];
	
	for( Box::pointIter mi = mStart; mi != mEnd; ++mi ) {
	  int m = *mi;
	  
	  //from m to n
	  Vec3 r0 = pn - sourcefield[m];
	  leafLevel.addCloseTransfer(n, m, r0, k_out);
	  
	  //from n to m
	  Vec3 r1 = -r0;
	  leafLevel.addCloseTransfer(m, n, r1, k_out);
	  
	}
      }
    }
    
    if( verbose >= 2 )
      cout << "   Close transfer functions Defined" << endl;
    
    if( verbose >= 2 )
      cout << "Initialization Done" << endl;
    
  }
  


  // Destructor
  MLFMM_Env::~MLFMM_Env() {
    for( int i = 0; i < (int) leaf_Translation.size(); i++ )
      delete leaf_Translation[i];
  }
  
  
  // Computes   omega = A*psi
  // void execute( const complex outgoing[], complex transfered[], complex kappa, Type type)
  void  MLFMM_Env::execute( const std::vector< std::vector<complex> >& PWE, std::vector< std::vector<complex> >& T_PWE, complex kappa, Type type)
  {
    int nLevels = tree.numLevels();
    if( Interface_Level != -1 ) getLevel(Interface_Level).initializeFields( type );
    Level& leafLevel = getLevel( nLevels );
    leafLevel.zeroFields();

    for( int i = 0; i < NScat; i++ )
      T_PWE[i].resize( global_index.size(), 0.);
    
    // ***** Bottom level translations *****
    double begin = clock();
    // For each box at leaf level
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;
      const Box::pointIter piEnd = b->pointIndex.end();
      
      // For each particle in box
      for( Box::pointIter pi = b->pointIndex.begin(); pi != piEnd; ++pi ) {
	int idx = *pi;

	// Proceed to translation
	leaf_Translation[idx]->TranslateUp(b->getMultipole(), PWE[idx], FORWARD);
	
      }
    }
    
    if( verbose >= 1 ) cout << "Lowest Level Multipoles Done : " <<  (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    



    
    // ***** Upward Pass *****

    begin = clock();
    for( int L = nLevels; L >= 3; --L ) 
      {
	Level& level = getLevel(L);
	Level& pLevel = getLevel(L-1);
	pLevel.zeroFields();
	
	// Translation up (for all boxes at level L)
	for( BoxIter bi = level.boxbegin(); bi != level.boxend(); ++bi ) {
	  Box* b = *bi;

	  if( type == FORWARD )
	    pLevel.getTranslationUp(b).TranslateUp(b, level.getInterp(), FORWARD);
	  else if( type == ADJOINT )
	    pLevel.getTranslationDown(b).TranslateDown(b, level.getInterp(), ADJOINT);
	}

      }
    if( verbose >= 1 )
      cout << "Upward Pass Done : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    


    // ***** Transfer pass - transfer expansions between far-neighbors *****
    begin = clock();
    for( int L = nLevels; L >= 2; --L ) 
      getLevel(L).applyTransferFunctions( type );      

    if( verbose >= 1 ) cout << "Transfer List Applied : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    
   

    // ***** Downward pass *****
    begin = clock();
    for( int L = 3; L <= nLevels; ++L ) 
      {

	Level& level = getLevel(L);
	Level& pLevel = getLevel(L-1);
		
	//Translation (for all boxes at level L)
	for( BoxIter bi = level.boxbegin(); bi != level.boxend(); ++bi ) {
	  Box* b = *bi;

	  if( type == FORWARD )
	    pLevel.getTranslationDown(b).TranslateDown(b, pLevel.getAnterp(), FORWARD);
	  else if( type == ADJOINT )
	    pLevel.getTranslationUp(b).TranslateUp(b, pLevel.getAnterp(), ADJOINT);
       
	}

      }


    if( verbose >= 1 ) cout << "Downward Pass Done : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    


    
    // ***** Translation down to observation points *****
    begin = clock();
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;
      
      // for all interior points
      const Box::pointIter nEnd = b->pointIndex.end();
      for( Box::pointIter ni = b->pointIndex.begin(); ni != nEnd; ++ni ) {
	int n = *ni;

	// Proceed to translation
	if( type == FORWARD )
	  leaf_Translation[n]->TranslateDown(T_PWE[n], b->getLocal(), FORWARD);
	else if( type == ADJOINT )
	  leaf_Translation[n]->TranslateUp(T_PWE[n], b->getMultipole(), ADJOINT);
	
      }
    }
    
    if( verbose >= 1 ) cout << "Local Integration Done : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;




    // ***** Close-term interactions *****
    begin = clock();
    leafLevel.applyClose( PWE, T_PWE, kappa, type );
    
    if( verbose >= 1 ) cout << "Close List Applied : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
  }







  // Main algorithm to determine transfer pairs and close pairs
  // Two boxes are transfers if their parents are neighbors and they aren't
  // Two boxes are close if they are leaf level and are neighbors  
  void MLFMM_Env::defineTransferAndClose( Box* b, int L )
  {
    // For all the children of b
    for( int k1 = 0; k1 < tree.BRANCH; ++k1 ) {
      // Get the k1th child of b
      Box* child1 = getBox( tree.child(b->n, k1), L+1 );
      if( child1 == NULL ) continue;

      // For all the other children of b
      for( int k2 = k1+1; k2 < tree.BRANCH; ++k2 ) {
	// Get the k2th child of b
	Box* child2 = getBox( tree.child(b->n, k2), L+1 );
	if( child2 == NULL ) continue;

	// Call the neighbor routine on each pair of children
	defineTransferAndClose( child1, child2, L+1 );
      }

      // If the child is not at the leaf level, recurse on the child
      if( L+1 != tree.numLevels() )
	defineTransferAndClose( child1, L+1 );
    }
  }

  // Secondary algorithm to determine transfer pairs and close pairs
  void MLFMM_Env::defineTransferAndClose( Box* b1, Box* b2, int L )
  {
    // If leaf level, these boxes are a close pair
    if( L == tree.numLevels() ) {
      getLevel(L).addClose( b1, b2 );
      return;
    }

    // For all the children of b1
    for( int k1 = 0; k1 < tree.BRANCH; ++k1 ) {
      // Get the k1th child of b1
      Box* child1 = getBox( tree.child(b1->n, k1), L+1 );
      if( child1 == NULL ) continue;

      // For all the children of b2
      for( int k2 = 0; k2 < tree.BRANCH; ++k2 ) {
	// Get the k2th child of b2
	Box* child2 = getBox( tree.child(b2->n, k2), L+1 );
	if( child2 == NULL ) continue;

	// Determine the transfer std::vector between the two children
	Vec3 r0 = child1->center - child2->center;
	if( r0.mag() > 1.8 * tree.getBoxSize(L+1) ) {
	  // These two are not neighbors so they are a transfer pair
	  getLevel(L+1).addTransfer( child1, child2 );
	} else {
	  // These two are neighbors so recurse
	  defineTransferAndClose( child1, child2, L+1 );
	}
      }
    }
  }


}

#endif
