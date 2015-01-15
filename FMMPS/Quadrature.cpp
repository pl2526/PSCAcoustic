/*
 *  Quadrature.cpp
 *  PSCAcoustic
 *
 *  Quadrature for High-frequency portion of code
 *
 *
 *  Copyright (C) 2014 Pierre-David Letourneau
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef QUADRATURE_FMM_CPP
#define QUADRATURE_FMM_CPP

#include "Quadrature.h"

namespace FMM {

  //---- Quadrature points ----//

  // Constructor  
  QuadraturePoint::QuadraturePoint( double theta_, double phi_, double w_ ) : theta(theta_), phi(phi_), w(w_), s(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) {}

  
  //---- Quadrature rows ----//

  // Constructor
  QuadratureRow::QuadratureRow( double theta_, int N_phi, double w0 ) 
    : theta(theta_), point(N_phi)
  {
    // Construct the quadrature points
    for( int n = 0; n < N_phi; ++n ) {
      double phi = n * 2*PI/N_phi;
      double w = w0 * 2*PI/N_phi;
      point[n] = QuadraturePoint(theta, phi, w); 
    }
  }

    
    //---- Quadrature Helmholtz ----//
    
    // Constructor
    //template <typename FUNC>
    QuadratureHelm::QuadratureHelm( HelmKernel& K ) : kappa(std::abs(K.kappa)) {}
  


    //---- Quadrature Uniform ----//
    
    // A uniform quadrature in phi and theta
    
  Quadrature_Uniform::Quadrature_Uniform(HelmKernel& K, double boxSize, double eps, complex k_, 
					 complex k_out_) : QuadratureHelm(K), k(k_), k_out(k_out_)
  {
    //double A = 1;
    double r0 = 2.5 * boxSize;
    double r = 1.5 * sqrt(3.)/2. * boxSize;
    double r_g = 1.5 * sqrt(3.)/2. * boxSize;
    double r_h = 1.5 * sqrt(2.)/2. * boxSize;
    
    /* Determine the Gegenbauer truncation (L) */
    L = get_Truncature(std::abs(k_out), r, r0, eps, 3) + 30; 
    //cout << "Gegenbauer series truncation : " <<  L << endl;
    int maxL = L;  
    
     // *** Compute N_theta ***
    int N_theta = get_N_theta(L, std::abs(k_out), r0, r_g,  boxSize, eps);
    cout << "N_theta : " << N_theta << endl;
    int N_rows = N_theta/2 + 1;


    
 
    //Vec3 r0(norm_r0, 0, 0);
    //std::vector< std::vector<complex> > TL = mod_transfer_fn(L, kappa, N_rows, N_phiT, r0);
    

    // *** Compute N_phi(theta) ***
    row.resize( N_rows );
    int N_phi;
    int N_phiT = (2*L + 2 + 2*LMAX);
    Vec3 vec_r0(r0, 0, 0);
    std::vector< std::vector<complex> > TL = mod_transfer_fn(L, k_out, N_rows, N_phiT, vec_r0);
    int size = 0;
    for( int n = 0; n < N_rows; ++n ) {
      double theta = n * PI/(N_rows-1);
      
      if( n == 0 || n == N_rows-1 ) {   // The Poles
	
	double w0 = 2*PI/N_theta;
	
	row[n] = QuadratureRow( theta, 1, w0 );
	size++;
      } else if( n <= N_rows/2 ) {      // Positive z row
	
	// The weight of this row, 2x for symmetry
	double w0 = 2 * 2*PI/N_theta;	  
	
	//double error;
	N_phi = get_N_phi(TL, n, N_phiT, maxL, theta, N_theta, std::abs(k_out), r_h, boxSize, eps);
	cout << "N_phi : "<< N_phi << endl;
	
	row[n] = QuadratureRow( theta, N_phi, w0 );	
	size += N_phi;  
      } else {                         // Negative z row
	
	// The weight of this row, 2x for symmetry
	double w0 = 2 * 2*PI/N_theta;
	int N_phi_local = row[N_rows-1-n].size();
	
	row[n] = QuadratureRow( theta, N_phi_local, w0 );
	size += N_phi;	  
      }
    }
   
    makePointerAccess();
    cout << "Quadrature size : " << size << endl;

  }

  
    

  
  
  //-----Auxiliary functions-----//
  
  double Quadrature_Uniform::g( int n, double C, double k_out, double r){
    double sum = 0;
    for( int k = -2*LMAX; k <= 2*LMAX; k++ )
      sum += std::abs( bessel_J(n-k, k_out*r) );    

    return C * sum;  
  }
  
  
  double Quadrature_Uniform::h( int n, double C, double k_out, double r, double theta){
    double val = 0;
    for( int i = 0; i <= std::max(0, n + 2*LMAX); i++ )
      val = std::max(val, std::abs( bessel_J(i, k_out*r*sin(theta)) ) );
	    
    return C * val;
  }
  


    
  // TODO: NEED NOT DEPEND ON GLOBAL CONSTANT
  int Quadrature_Uniform::get_N_theta( int L, double k_out, double r0, double r,
				       double boxSize, double eps )
  {
    int N_rows_max = L;
    int N_theta_max = 8*N_rows_max+1;


    // *** Get the absolute value of the Fourier coefficients of the transfer function ***

    // Get the transfer matrix for phi = 0, theta = [0,2PI)
    Vec3 vec_r0(0,0,r0);
    std::vector< std::vector<complex> > T = mod_transfer_fn(L, k_out, N_rows_max, 2, vec_r0);
    
    // Copy real-space values into a single std::vector t(theta)
    std::vector<complex> t(N_theta_max);
    t[0] = T[0][0];                        // North Pole
    for( int k = 1; k < N_rows_max; ++k ) {
      t[k] = T[k][0];
      t[N_theta_max-k] = T[k][1];
    }
    
    // Fourier coefficients
    ifft( &t[0], N_theta_max );
    std::vector<double> absT(N_theta_max);
    for( int k = 0; k < N_theta_max; ++k ) 
      absT[k] = std::abs(t[k]);
    
    

    // *** Normalizing constant ***
    double C = 0;
    for( int l = 0; l <= LMAX; l++){
      for( int l_p = 0; l_p <= LMAX; l_p++){
       	
	double x =  1./(2*PI*PI) * std::abs(Scatterer::Mie(K, k_out, RADIUS, RHO, 0, l)); 
	x /= std::abs(Amos_sf_hankel_1(l_p, 1.1 * k_out*RADIUS));

	C = std::max(C, x);
      }
    }


    // *** Search for N_theta ***
    //for( int N_theta = 2*L; N_theta < N_theta_max; N_theta += 2 ) {
    for( int N_theta = 2; N_theta < N_theta_max; N_theta += 2 ) {      

      // Get the truncation error sum_{n >= N_theta/2} + sum_{n <= -N_theta/2}
      double error_trunc = 0;
      for( int n = N_theta/2; n < N_theta_max; ++n )
	error_trunc += 2 * absT[n] * g(n, C, std::abs(k_out), r);
      
      // Get the aliasing error (Only compute p=1)
      double error_alias = 2 * absT[0] * g(N_theta, C, std::abs(k_out), r);
      for( int n = 1; n < N_theta/2; ++n )
	error_alias += 2 * absT[n] * g(n + N_theta, C, std::abs(k_out), r);
      
      double error = 4*PI*PI * 4*PI/std::abs(k_out) * ( error_alias + error_trunc );
      if( error < eps ) return N_theta;
    }
    
    cerr << "Error: N_theta convergence failure" << endl;
    return N_theta_max;
  }
  
  
  int Quadrature_Uniform::get_N_phi( std::vector< std::vector<complex> >& TL, 
				     int n,  int N_phiT, double maxL, double theta, 
				     int N_theta, double k_out, double r, double boxSize, double eps)
  {
    int N_phi;
    
    // *** Get the absolute value of the Fourier coefficients of the transfer function ***
    ifft( &TL[n][0], N_phiT );	
    
    std::vector<double> absT(L+1);
    for( int m = 0; m <= L; ++m ) 
      absT[m] = std::abs( TL[n][m] );
    
    


    // *** Normalizing constant ***
    double C = 0;
    for( int l = 0; l <= LMAX; l++){
      for( int l_p = 0; l_p <= LMAX; l_p++){
	
	// Normalize
	double x = sqrt(l_p+1.)*sqrt(l+1.)/(4.*PI) * std::abs(Scatterer::Mie(K, k_out, RADIUS, RHO, 0, l)); 
	x /= std::abs(Amos_sf_hankel_1(l_p, 1.1 * k_out*RADIUS));
	
	C = std::max(C, x);
      }
    }




    // *** Compute N_phi by searching for error ***
    for( N_phi = 4; N_phi <= ceil4( 2*maxL+1+2*LMAX ); N_phi += 4 ) {
      
      // Get the truncation error sum_{|k| >= N_phi/2}
      double error_trunc = 0;
      for( int k = N_phi/2; k <= L; ++k )
	error_trunc += 2 * absT[k] * h(k, C, std::abs(k_out), r, theta);
      
      // Get the aliasing error
      double error_alias = 2 * absT[0] * h(N_phi, C, std::abs(k_out), r, theta);
      for( int k = 1; k <= N_phi/2-1; ++k ) {
	error_alias += 2 * absT[k] * h(N_phi-k, C, std::abs(k_out), r, theta);
      }
      
      double error = 4*PI*PI * 4*PI/std::abs(k_out) * ( error_alias + error_trunc ); 
      
      if( error < eps ) break;
    }
    
    return N_phi;
  }
    
  
}

  
#endif
