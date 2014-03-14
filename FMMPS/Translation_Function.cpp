#ifndef TRANSLATION_FUNCTION_FMM_CPP
#define TRANSLATION_FUNCTION_FMM_CPP

#include "Translation_Function.h"

namespace FMM {


  // Initialize static index
  Indexing Translation_Function::global_index(LMAX);
  
  // TODO: Still needed??
  
  //Translate partial wave expansion of each scatterer to center of its box at the leaf level
  //from vec2 to vec1
  /* void Translation_Function::ManualTranslation_Direct( const Vec3& r, complex k_out,
						       std::vector<complex>& vec1,  Indexing* index1, 
						       const std::vector<complex>& vec2, Indexing* index2, bool reg) 
  {
    assert( (int) vec1.size() == (int) index1->size() );
    assert( (int) vec2.size() == (int) index2->size() );
    
    //TODO : merge Pvec and Vec3
    double rad, theta, phi;
    Cart2Sf(r.x, r.y, r.z, rad, theta, phi);
    
    
    PointAndShoot PS(rad, k_out, theta, phi, index1, index2, reg);
    vec1 = PS.Apply(vec2, FORWARD);
  }
  


    //Translate Far-field signature from vec2 to vec1
  void Translation_Function::ManualTranslation_HF(const Vec3& r, 
						  std::vector<complex>& vec1, const std::vector<complex>& vec2, 
						  Quadrature* quad) 
  {
    assert( (int) vec1.size() == (int) quad->size() );
    assert( (int) vec2.size() == (int) quad->size() );
    
    double kappa = quad->kappa;
    int K = quad->size();
    for( int k = 0; k < K; ++k ) {
      const double kr = kappa * r.dot( quad->getPoint(k).s );
      vec1[k] = vec2[k] * complex( cos(kr), sin(kr) );
    }
  }

  */


  
  //-- Translation in the high-frequency regime --//
 
  
  // Precomputes the Translation function for std::vector r
  HF_Translation_Function::HF_Translation_Function(const Vec3& r, complex k_out, bool INTF_, 
						   SH_Transform* SH_T_, Quadrature* quad) 
    : Translation_Function(), INTF(INTF_), SH_T(SH_T_)
  {
    int K = quad->size();
    C.resize(K);
    
    for( int k = 0; k < K; ++k ) {
      const Vec3& s = quad->getPoint(k).s;
      C[k] = exp( CI * k_out * r.dot(s) );
    }
    
  }
  
  
  // ***Note : Reterpolators always go from q1 to q2
    // Translate multipole expansion up from child (M) to parent (pM)
  void HF_Translation_Function::TranslateUp(Box* b, Interpolator* interp, Type type)
  {
    // From level to pLevel
    level_quad = interp->q1;
    pLevel_quad = interp->q2;
    
    // Interpolate to the parent level L-1
    std::vector<complex> M;
    if( INTF )
      M = SH_T->Forward( b->getMultipole() );
    else
      M = b->getMultipole();
    
    
    assert( (int) M.size() == level_quad->size() );
    pM.assign(pLevel_quad->size(), 0.);
    interp->apply( M, level_quad, pM, pLevel_quad );
    
    // Multiply by the transfer function and accumulate into parent
    int N = b->parent->getMultipole().size();
    assert( N == (int) C.size() );
    assert( pM.size() == C.size() );
    for( int i = 0; i < N; i++)
      b->parent->getMultipole()[i] += C[i] * pM[i];
  }
  

  
  // Translate multipole expansion down from parent (pM) to child (M)
  void  HF_Translation_Function::TranslateDown(Box* b, Anterpolator* anterp, Type type)
  {
    // From pLevel to level
    pLevel_quad = anterp->q1;
    level_quad = anterp->q2;
    
    M.assign( level_quad->size(), 0.);
    pM.assign( pLevel_quad->size(), 0.);
    
    // Multiply parent's local expansion by translation function
    int N = b->parent->getLocal().size();
    assert( N == (int) pM.size() );
    assert( N == (int) C.size() );
    for( int i = 0; i < N; i++ )
      pM[i] = C[i] * b->parent->getLocal()[i];
    
    // Anterpolate
    anterp->apply( pM, pLevel_quad, M, level_quad);
    
    // If interface level, proceed to SHT
    if( INTF )
      M = SH_T->Backward( M );
    
    // Accumulate into the box's local expansion
    int K = b->getLocal().size(); 
    assert( K == (int) M.size() );
    for( int i = 0; i < K; i++)
      b->getLocal()[i] += M[i];	
  }
  
  
  
  
  
  
  
  
  //-- Translation in the low-frequency(LF) regime--//
  
  // Constructor
  Direct_Translation_Function::Direct_Translation_Function( Vec3 r_, complex k_out, Indexing* pIndex, 
							    Indexing* index, Indexing* cIndex)
    : Translation_Function(), r(r_)
  {
    double rad, theta, phi;
    Cart2Sf(r.x, r.y, r.z, rad, theta, phi);
    
    // From child to current
    PS_Translation_Up = new PointAndShoot(rad, k_out, theta, phi, index, cIndex, true);
    // From current to child
    PS_Translation_Down = new PointAndShoot(rad, k_out, theta, phi, cIndex, index, true);
  }
  
  // Destructor
  Direct_Translation_Function::~Direct_Translation_Function()
  {
    delete PS_Translation_Up;
    delete PS_Translation_Down;
  }
  
  
  void Direct_Translation_Function::TranslateUp(Box* b,  Interpolator* interp, Type type)
  {
    // Perform translation
    std::vector<complex> pM = PS_Translation_Up->Apply(b->getMultipole(), FORWARD);
    
    // Accumulate into parent expansion
    assert( pM.size() == b->parent->getMultipole().size() );
    for( int i = 0; i < (int) pM.size(); i++ )
      b->parent->getMultipole()[i] += pM[i];
  }
  

  void Direct_Translation_Function::TranslateDown(Box* b, Anterpolator* anterp, Type type)
  {
    // Perform translation
    std::vector<complex> M;
    M = PS_Translation_Down->Apply(b->parent->getLocal() , FORWARD );
    
    // Accumulate into child expansion
    assert( M.size() == b->getLocal().size() );
    for( int i = 0; i < (int) M.size(); i++ )
      b->getLocal()[i] += M[i]; 
  }
  
  
  
  
  
  
  
  //-- Translation at the leaf level (last level of octree) --//
  // TODO : Improve structure
  
  // Precomputes the Translation function for std::vector r
  Leaf_HF_Translation_Function::Leaf_HF_Translation_Function(const Vec3& r, complex k_out, Quadrature* leaf_quad) 
    : Leaf_Translation_Function()
  {
    // TODO: delete in destructor
    // Construct SHT
    SH_T = new SH_Transform( leaf_quad, &Translation_Function::global_index );

    // Compute translation elements
    C.resize( leaf_quad->size() );
    for( int k = 0; k < (int) leaf_quad->size(); ++k ) {
      const Vec3& s = leaf_quad->getPoint(k).s;
      C[k] = exp( CI * k_out * r.dot(s) );
    }
    
  }
  
  
  void Leaf_HF_Translation_Function::TranslateUp(std::vector<complex>& To, std::vector<complex> From, Type type)
  {
    // Perform Forward SHT
    From = SH_T->Forward( From );
    
    // Multiply by the transfer function and accumulate into M
    int K = C.size();
    assert( (int) To.size() == K );
    for( int i = 0; i < K; i++)
      To[i] += C[i] * From[i];
  }
  
  
  void Leaf_HF_Translation_Function::TranslateDown(std::vector<complex>& To, std::vector<complex> From, Type type)
  {
    //int K = C.size(); 
    assert( C.size() == From.size() );
    for( int i = 0; i < (int) C.size(); i++)
      From[i] = conj(C[i]) * From[i];
    
    // Perform Backward SHT
    From = SH_T->Backward( From );
    
    assert( From.size() == To.size() );	
    for( int i = 0; i < (int) To.size(); i++)
      To[i] += From[i];
  }
  
  
  
  // Low-frequency regime
 
  // Constructor
  Leaf_Direct_Translation_Function::Leaf_Direct_Translation_Function( Vec3 r, complex k_out, Indexing* leafIndex ) 
    : Leaf_Translation_Function()
  {
    double rad, theta, phi;
    Cart2Sf(r.x, r.y, r.z, rad, theta, phi);
    T_up = new PointAndShoot(rad, k_out, theta, phi, leafIndex, &Translation_Function::global_index, true);
    
    Cart2Sf(-r.x, -r.y, -r.z, rad, theta, phi);
    T_down = new PointAndShoot(rad, k_out, theta, phi, &Translation_Function::global_index, leafIndex, true);
  }
  
  
  // Note: leafQuad and leafIndex are freed by the leaf level itself
  // Destructor
  Leaf_Direct_Translation_Function::~Leaf_Direct_Translation_Function()
  {
    delete T_up;
    delete T_down;
  }
  
  
  void Leaf_Direct_Translation_Function::TranslateUp(std::vector<complex>& To, std::vector<complex> From, Type type)
  {
    // Perform translation up
    vec =  T_up->Apply(From, type );
    assert( vec.size() == To.size() );
    for( int i = 0; i < (int) To.size(); i++)
      To[i] += vec[i];
  }
  
  
  void Leaf_Direct_Translation_Function::TranslateDown(std::vector<complex>& To, std::vector<complex> From, Type type)
  {
    // Perform translation down
    vec =  T_down->Apply(From, type );
    assert( vec.size() == To.size() );
    for( int i = 0; i < (int) To.size(); i++)
      To[i] += vec[i];
  }
  
  
   
}

#endif