#ifndef TRANSFER_FUNCTION_FMM_H
#define TRANSFER_FUNCTION_FMM_H

#include "Quadrature.h"
#include "TransferUtils.h"

namespace FMM {
  

  //---- Super class for transfer function ----//
  class Transfer_Function
  {
    
  public:

    static Indexing global_index;

    // Constructor/Destructor
    Transfer_Function(){}
    virtual ~Transfer_Function(){}
    
    // Execution
    // Transfer from vec1 to vec2
    virtual void Transfer( std::vector<complex>& vec1, std::vector<complex>& vec2, Type type ) { return; }  
  };
  
  
  
  //---- High-frequency transfer function ----//
  class HF_Transfer_Function : public Transfer_Function
  {
    std::vector<complex> C;
    complex k_out;

  public:
    
    // Constructor
    HF_Transfer_Function(Quadrature* quad, const Vec3& r0, complex k_out_);
    
    //Destructor
    ~HF_Transfer_Function(){}

    // Execution
    // from vecFrom to vecTo
    virtual void Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type );
  };




  //---- Low-frequency transfer function ----//
  class Direct_Transfer_Function : public Transfer_Function 
  {
 
    PointAndShoot* PS_Transfer;

  public:
    // Constructor/Destructor
    //Direct_Transfer_Function( Indexing* index, Vec3 r, complex k_out);
    Direct_Transfer_Function( Indexing* index, Vec3 r, complex k_out);
    ~Direct_Transfer_Function();

    // From vecFrom to vecTo
    virtual void Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type );
  };

}

#endif