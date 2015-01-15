#ifndef LEVEL_FMM_H
#define LEVEL_FMM_H

#include "General.h"
#include "Transfer_Function.h"
#include "Transfer_Vector.h"
#include "Translation_Function.h"
#include "NFunction.h"
#include "TransferUtils.h"
#include "Quadrature.h"
#include "../L_choice.h"


namespace FMM {

  
  typedef list<Trans_Vector> TVList;
  typedef TVList::iterator TVIter;
  
  typedef list<Close_Vector> CVList;
  typedef CVList::iterator CVIter;
  
  typedef std::vector<Transfer_Function*> TFList;
  typedef TFList::iterator TFIter;
  
  typedef std::vector<Translation_Function*> TLList;
  typedef TLList::iterator TLIter;
  
  typedef list<Box*> BoxSet;
  typedef BoxSet::iterator BoxIter;
  
  //template <int DIM>
    class Level
    {  
    public:
      const static int BRANCH = 1 << DIM;   // The branching factor = 2^DIM

      bool HF;         // Level belongs to HF regime or not
      bool INTF;       // whether level lies at interface between HF and LF regime
      bool LEAF;       // whether level is leaf level

      int L;                   // Max order for spherical wave functions
      Indexing* index;         //Indexing for spherical wave functions up to order L
      static Indexing global_index;
      
      double boxSize;          // Size of a box this levels contains
      complex k_out;           // Background wave number
      complex k;               // Scatterer wave number
      
      BoxSet boxList;          // The list of nonempty boxes this level contains
      
      Quadrature* quad;        // Quadrature for numerical functions of this level
      NFunction scratchSpace;  // Scratch space of size quadrature for computation
      
      TVList transList;        // The list of transfer std::vectors
      CVList closeList;        // The list of close std::vectors (necessary?)
      
      TFList transferFunc;     // The list of transfer functions (necessary?)
      TLList translationFunc;  // The list of translation functions

      std::vector< std::vector<int> > closeIdx;  // The index of the particles associated w/ a close transfer
      std::vector<Vec3> closeVec;           // The transfer std::vectors
      TFList closeFunc;        // The list of close transfers (only at leaf level)
      
      Interpolator* upReterp;
      Anterpolator* downReterp;

      SH_Transform* SH_T;
      
    public:
      
      // Constructors
      Level( double boxSize_ = 0, complex k_out_ = K_OUT, complex k_ = K, double cutoff = 2*PI*OMEGA );

      // Destructor
      ~Level();


      //-- Auxilliary functions --//      
      
      inline std::vector< std::vector<int> > getCloseIdx() { return closeIdx; }
      inline double getBoxSize(){ return boxSize; }
      inline void addBox( Box* b ){ boxList.push_back(b); }
      inline Indexing* getIndex(){ return index; }
      inline Quadrature* getQuad(){ return quad; }
      inline int getIdxSize(){ return index->size(); }

      /* inline NFunction* getScratch()
      {
	return &scratchSpace;
	}*/

      void setL();
      void initializeFields( Type type = FORWARD );
      void zeroFields();




      //--- Quadrature ---//

      void defineQuadrature( HelmKernel KERNEL, double eps);


      //--- Spherical Harmonics Transform ---//

      // Getter
      inline SH_Transform* getSH_T(){ return SH_T; }

      // Create instance of Spherical Harmonics Transform
      inline void defineSH_T(Quadrature* quad_lower, Indexing* index_lower){
	SH_T = new SH_Transform(quad_lower, index_lower);	
      }
      
      // Compute actual transforms
      void computeSH_T_Forward( Type type );
      void computeSH_T_Backward( Type type );




      //--- Translation and Transfer ---//

      void defineTransfers( double eps, complex k_out );
      void addTransfer( Box* b1, Box* b2 );
      void addClose( Box* b1, Box* b2 );

      // Define close transfers (only at leaf level)
      void addCloseTransfer( int& idxTo, int& idxFrom, Vec3& r0, complex k_out );
      
      // Define the translation functions from a lower level
      void defineTranslations(Indexing* pIndex, Indexing* cIndex, complex k_out);


      // Get the translation function from a child box b
      // to its parent box in this level
      inline Translation_Function& getTranslationUp( Box* b )
      {
	return *translationFunc[b->n & (BRANCH-1)];
      }
      
      // Get the translation function from a parent box to a child box b
      // in this level
      inline Translation_Function& getTranslationDown( Box* b )
      {
	return *translationFunc[(b->n & (BRANCH-1)) ^ (BRANCH-1)];
      }
      



      //--- Interpolation ---//

      // TODO :Should follow same convention as translation functions
      // Interpolator: interpolation from current level to upper/parent level
      // Anterpolator: interpolation from current level to lower/child level
      inline void defineInterp( Quadrature* quad_Upper)
      {
	upReterp = new Interpolator( quad, quad_Upper);
      }
      
      inline Interpolator* getInterp()
      {
	return upReterp;
      }
      
      
      // Define an anterpolator from the quadrature of this level
      // to the quadrature of level qb
      inline void defineAnterp( Quadrature* quad_Lower )
      {
	downReterp = new Anterpolator( quad, quad_Lower);
      }
      
      inline Anterpolator* getAnterp()
      {
	return downReterp;
      }



      //--- Execution ---//
      
      void applyTransferFunctions( Type type );
      void applyClose(const std::vector< std::vector<complex> >& PWE, 
		      std::vector< std::vector<complex> >& T_PWE, 
		      complex k_out, Type type );
      
      // Get box indices
      inline BoxIter boxbegin(){ return boxList.begin(); }
      inline BoxIter boxend(){ return boxList.end(); }

    };
  
}
#endif
