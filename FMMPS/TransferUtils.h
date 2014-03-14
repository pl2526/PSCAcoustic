/*
 *  WignerIntegrals.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 4/25/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Computation/Storage of integrals of the form \int Y_u^v Y_m^n Y_p^q (Wigner-like integrals) appearing in the     *  transfer operator for spherical wave functions. 
 */
#ifndef TU_FMM_H
#define TU_FMM_H

#include "General.h"
#include "../General.h"
#include "../Indexing.h"

namespace FMM {

// TODO: Should not have different instances of Index
//Indexing Index(LMAX);

//**** Point and Shoot ****//


 
//--- Rotation about z-axis ---//

 class Z_Rotation
 {
   std::vector<complex> matrix;
 public:
   
   enum Sign{POS, NEG};
   
   // Constructor/Destructor
   Z_Rotation( Indexing* index, double theta ); 
   ~Z_Rotation(){}
   
   
 private:
   
   static void Compute( std::vector<complex>& matrix, Indexing* index, double theta );
   
 public:
   // Proceed to translation in-place
   void Apply( std::vector<complex>& vec, Indexing* index, Sign S = POS  );   
 };







  //--- S Matrix ---//

  struct S_Rotation 
  {
    static complex **S_matrix;
    static Indexing* index;
    
    // Constructor/Destructor
    S_Rotation(){}
    ~S_Rotation(){}
 
    
    // Initialization/Construction of the rotation matrix elements
    static void S_Init( int L );

    // Needed to delete static pointer
    static void S_Delete();

    // Some auxilliary inline functions    
    static inline double a(int l, int m, int m_p){ return sqrt( (double) (l+m) * (l-m) / ((l+m_p) * (l-m_p)) ); }
    static inline double b(int l, int m, int m_p){ return sqrt( (double) (l+m) * (l+m-1) / (2*(l+m_p) * (l-m_p)) ); }
    static inline double c(int l, int m, int m_p){ return sqrt( (double) 2*(l+m) * (l-m) / ((l+m_p) * (l+m_p-1)) ); }
    static inline double d(int l, int m, int m_p){ return sqrt( (double) (l+m) * (l+m-1) / ((l+m_p) * (l+m_p-1)) ); }
    
       
  public:

    // Application operator
    void Apply( std::vector<complex>& vec, int L, Type type = FORWARD);

  };

  
  
  
  
  
  
  //--- Transfer along z-axis ---//
  
  //Transfer from 2 to 1 along z-axis
  class Z_transfer
  {
    Indexing* index1;
    Indexing* index2;
    int L_1, L_2, L;
    static Indexing* index;
    
    std::vector<complex> coeff;
    
  public:
    
    // Constructor/Destructor
    Z_transfer( double R, complex k_out, Indexing* index1_, Indexing* index2_, bool reg);
    ~Z_transfer(){}
    
    // Initialization/Destruction of Z_transfer matrix
    static void Z_Init( int L );
    static void Z_Delete();
    
    
    // Recursive computation of transfer  and translation coefficients
    // For details, see Gumerov, Duraiswami : "Fast, Exact, and Stable Computation of Multipole Translation
    // and Rotation Coefficients for the 3-D Helmholtz equation"
    // Matlab Prototype : GumerovTransferCoefficients.m
    
    // Goes from (l,m) to (l',m')
    void Compute(SpMatrix& coeff2, int L, double r, complex k, bool reg);

    void Apply( std::vector<complex>& vec, Type type = FORWARD );
  };


  
  
  
  
  // Accelerated Translation/Transfer from 2 to 1
  class PointAndShoot 
  {
    // Z-rotation matrices
    Z_Rotation* Z_phi;      // Phi-rotation along z-axis
    Z_Rotation* Z_theta;    // Theta-rotation along z-axis
    
    // TODO : should be static and have enough elements to accomodate any relevant index
    // S-rotation matrices
    S_Rotation* S_F ; // Forward S rotation matrix
    S_Rotation* S_B; // Backward S rotation matrix
    //S_Rotation S;
    
    // Z-transfer along z-axis
    Z_transfer* Z_T; // Transfer along the z-axis
    
    Indexing* index1;
    Indexing* index2;
    Indexing* index; 
    
  public:
    
    // Constructor/Destructor
    // (Accelerated Translation/Transfer from 2 to 1)
    PointAndShoot(double r, complex k_out, double theta, double phi, Indexing* index1_, Indexing* index2_, bool reg);
    ~PointAndShoot();
    
    // Application of operator
    // Note : pass by value since rotations act in place
    std::vector<complex> Apply( std::vector<complex>& vec, Type type = FORWARD );
  };

  
  
  
  // Computes full transfer matrix from 2 to 1 along the std::vector (r,theta,phi)   
  // Recursive computation of transfer  and translation coefficients
  // For details, see Gumerov, Duraiswami : "Fast, Exact, and Stable Computation of Multipole Translation
  // and Rotation Coefficients for the 3-D Helmholtz equation"
  // Matlab Prototype : GumerovTransferCoefficients.m
  class Gumerov_TransCoeff
  {
    
  public:
    // TODO: use Eigen
    arma::Mat<complex> TCMatrix;
    int L;
    
    // Constructor
    Gumerov_TransCoeff(double r, complex k_out, double theta, double phi, int L_, bool reg);
    
    // Compute entries (goes from (l,m) to (l',m'))
    static arma::Mat<complex> Compute(int L, double r, double theta, double phi, complex k, bool reg);
    
    // TODO: Needed?
    // Apply
    //std::vector<complex> Apply( std::vector<complex> vec, Type type = FORWARD );
    
    
    
    
    // Some inline functions
    static inline double a(int n,int m){
      double entry;
      ( n >= std::abs(m) ) ? entry = sqrt( (n+1.+std::abs(m)) * (n+1.-std::abs(m)) / ((2.*n+1.)*(2.*n+3.)) ) : entry = 0. ;
      return entry;
    }
    
    static inline double b(int n,int m){
      double entry;
      if( 0 <= m && m <= n ){
	entry = sqrt( (n-m-1.)*(n-m) / ((2.*n-1.)*(2.*n+1.)) );
      } else if( -n <= m && m < 0 ){
	entry = -sqrt( (n-m-1.)*(n-m) / ((2.*n-1.)*(2.*n+1.)) );
      } else {
	entry = 0;
      }
      
      return entry;
    }
    
    
  };


}


#endif
