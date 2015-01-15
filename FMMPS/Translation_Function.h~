#ifndef TRANSLATION_FUNCTION_FMM_H
#define TRANSLATION_FUNCTION_FMM_H

#include "Reterpolator.h"
#include "Quadrature.h"
#include "SphericalHarmonicsTransform.h" 
#include "Box.h"
#include "TransferUtils.h"

#include "../Constants.h"
#include "../Coordinates.h"
#include "../Indexing.h"


namespace FMM {

  // *** TREE CORE TRANSLATIONS ***

  // Super-class for translations
  class Translation_Function{

  public:
    static Indexing global_index;

    Translation_Function(){}
    virtual ~Translation_Function(){}
    
    
    //TODO : Are these manual translations still needed?
    //TODO : If so update implementation for adjoint
    
    //Translate partial wave expansion of each scatterer to center of its box at the leaf level
    //from vec2 to vec1
    void static ManualTranslation_Direct( const Vec3& r, complex k_out,
					  std::vector<complex>& vec1,  Indexing* index1, 
					  const std::vector<complex>& vec2, Indexing* index2, bool reg);

    //Translate Far-field signature from vec2 to vec1
    void static ManualTranslation_HF(const Vec3& r, 
				     std::vector<complex>& vec1, const std::vector<complex>& vec2, 
				     Quadrature* quad);
    
    // Prototypes for virtual functions
    virtual void TranslateUp(Box* b, Interpolator* interp = NULL, Type type = FORWARD) { return; }
    virtual void TranslateDown(Box* b,  Anterpolator* anterp = NULL, Type type = FORWARD) { return; }
  };
  



  
  //-- Translation in the high-frequency regime --//
  class HF_Translation_Function : public Translation_Function
  {
   
    std::vector<complex> C;

    Quadrature* level_quad;
    Quadrature* pLevel_quad;
    Quadrature* cLevel_quad;

    std::vector<complex> M;    // Scratch space for box expansion
    std::vector<complex> pM;   // Scratch space for box's parent expansion
    std::vector<complex> cM;   // Scract space for child's parent expansion

    bool INTF;
    SH_Transform* SH_T;
  public:
    
    // Constructor
    HF_Translation_Function(const Vec3& r, complex k_out, bool INTF_, SH_Transform* SH_T_, Quadrature* quad = NULL);

    // Destructor
    ~HF_Translation_Function(){}

    
    // *** Note : Reterpolators always go from q1 to q2 ***
    // Translate multipole expansion up from child (M) to parent (pM)
    virtual void TranslateUp(Box* b, Interpolator* interp = NULL, Type type = FORWARD);
    
    // Translate multipole expansion down from parent (pM) to child (M)
    virtual void TranslateDown(Box* b, Anterpolator* anterp = NULL, Type type = FORWARD);
  };
  






  //-- Translation in the low-frequency(LF) regime--//
  class Direct_Translation_Function : public Translation_Function {

    Vec3 r;
    PointAndShoot *PS_Translation_Up;
    PointAndShoot *PS_Translation_Down;

  public:    
    
    // Constructor
    Direct_Translation_Function( Vec3 r_, complex k_out, Indexing* pIndex, Indexing* index, Indexing* cIndex);

    // Destructor
    ~Direct_Translation_Function();

    virtual void TranslateUp(Box* b,  Interpolator* interp = NULL, Type type = FORWARD);
    virtual void TranslateDown(Box* b, Anterpolator* anterp = NULL, Type type = FORWARD);
  };





  // *** LEAF LEVEL TRANSLATIONS ***


  // TODO : Improve structure
  // ****Translation function at bottom of tree

  // Super-class for leaf translation
  class Leaf_Translation_Function
  {
    
  public:
    Leaf_Translation_Function(){}
    virtual ~Leaf_Translation_Function(){}
    
    virtual void TranslateUp(std::vector<complex>& M, std::vector<complex> expansion, Type type = FORWARD) { return; }
    virtual void TranslateDown(std::vector<complex>& M, std::vector<complex> expansion, Type type = FORWARD) { return; }
  };
  


  //-- Leaf translation in the high-frequency regime --//
  class Leaf_HF_Translation_Function : public Leaf_Translation_Function
  {

    std::vector<complex> C; // Storage of translation coefficients
    SH_Transform* SH_T;        // Spherical harmonics transform

  public:
    
    // Constructor
    Leaf_HF_Translation_Function(const Vec3& r, complex k_out, Quadrature* leaf_quad);

    // Destructor
    ~Leaf_HF_Translation_Function(){}

    virtual void TranslateUp(std::vector<complex>& To, std::vector<complex> From, Type type = FORWARD);
    virtual void TranslateDown(std::vector<complex>& To, std::vector<complex> From, Type type = FORWARD);
    
  };
  


  //-- Leaf translation in the low-frequency regime --//
  class Leaf_Direct_Translation_Function : public Leaf_Translation_Function
  {

    PointAndShoot *T_up;         // Translation up
    PointAndShoot *T_down;       // Translation down
    std::vector<complex> vec;

  public:    
    
    // Constructor
    Leaf_Direct_Translation_Function( Vec3 r, complex k_out, Indexing* leafIndex );

    // Note: leafQuad and leafIndex are freed by the leaf level itself
    // Destructor
    ~Leaf_Direct_Translation_Function();

    virtual void TranslateUp(std::vector<complex>& To, std::vector<complex> From, Type type = FORWARD);
    virtual void TranslateDown(std::vector<complex>& To, std::vector<complex> From, Type type = FORWARD);
  };



  
}

#endif
