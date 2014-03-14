#ifndef QUADRATURE_FMM_CPP
#define QUADRATURE_FMM_CPP

#include "Quadrature.h"

namespace FMM {

  //---- Quadrature points ----//

  // Constructor  
  QuadraturePoint::QuadraturePoint( double theta_, double phi_, double w_ ) 
  : theta(theta_), phi(phi_), w(w_), s(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) {}

  
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
    
    Quadrature_Uniform::Quadrature_Uniform(HelmKernel& K, double boxSize, double eps, complex k_, complex k_out_) 
      : QuadratureHelm(K), k(k_), k_out(k_out_)
    {

      double A = 1;
      
      // TODO: FIX
      //for( int l = 0; l <= LMAX; l++ )
      //    A = std::max(A, std::abs(Mie(k, k_out, RADIUS, RHO, l) * Amos_sf_hankel_1(l, k_out*RADIUS)) );
    
    
    double delta = 0.8;
    // Characteristic lengths
    double norm_r0 = 3. * boxSize;;//2 * boxSize;
    double norm_r = delta*sqrt(3) * boxSize;
    
    /* Determine the Gegenbauer truncation (L) */
    L = get_Truncature(kappa, norm_r, norm_r0, eps/(2*LMAX+1), 3)+5; 
    cout << "Gegenbauer series truncation : " <<  L << endl;
    int maxL = L;  
    
    /* Determine the number of quadrature rows (N_rows = N_theta/2 + 1) */
    int N_theta = get_N_theta(L, kappa, norm_r0, boxSize, A, eps);
    int N_rows = N_theta/2 + 1;
    int N_phiT = (2*L + 2 + 2*LMAX);
    
    cout << "N_theta : " << N_theta << endl;
    
    row.resize( N_rows );
    Vec3 r0(norm_r0, 0, 0);
    std::vector< std::vector<complex> > TL = mod_transfer_fn(L, kappa, N_rows, N_phiT, r0);
    
    // Prune N_phi(theta) if possible
    int N_phi;
    int size = 0;
    for( int n = 0; n < N_rows; ++n ) {
      double theta = n * PI/(N_rows-1);
      
      if( n == 0 || n == N_rows-1 ) {   // The Poles
	
	double w0 = 2*PI/N_theta;
	N_phi = 1;
	
	row[n] = QuadratureRow( theta, N_phi, w0 );
	size++;
      } else if( n <= N_rows/2 ) {      // Positive z row
	
	// The weight of this row, 2x for symmetry
	double w0 = 2 * 2*PI/N_theta;	  
	
	double error;
	N_phi = get_N_phi(error, TL, n, N_phiT, maxL, theta, N_theta, boxSize, A, eps);
	
	if( verbose == 2)
	  cout << "N_phi : "<< N_phi << endl;
	
	row[n] = QuadratureRow( theta, N_phi, w0 );	
	size += N_phi;  
      } else {                         // Negative z row
	
	// The weight of this row, 2x for symmetry
	double w0 = 2 * 2*PI/N_theta;
	N_phi = row[N_rows-1-n].size();
	
	row[n] = QuadratureRow( theta, N_phi, w0 );
	size += N_phi;	  
      }
    }
    
    makePointerAccess();
    cout << "Quadrature size : " << size << endl;
  }
  
    


  
  //-----Auxiliary functions-----//
  
  double Quadrature_Uniform::Gamma( int n, double kb, double A){
    double val = 0;
    for( int l_p = 0; l_p <= LMAX; l_p++){
      
      double sum = 0;
      for( int k = -(LMAX+l_p); k <= (LMAX+l_p); k++ )
	sum += std::abs( bessel_J(n-k, kb) );
      sum *= A * 1./(2*PI*PI) * 1./std::abs(Amos_sf_hankel_1(l_p, k_out*RADIUS));
      
      val = std::max(val, sum);
    }
    return val;  
  }
  
  
  double Quadrature_Uniform::Eta( int n, double kb, double theta, double A){
    
    double val = 0;
    for( int l_p = 0; l_p <= LMAX; l_p++){
      int idx = std::max(0, n-(LMAX+l_p));

      double C = sqrt(l_p+1)*sqrt(LMAX+1)/(4*PI) * 1./std::abs(Amos_sf_hankel_1(l_p, k_out*RADIUS));
      
      val = std::max(val, A*C*std::abs( bessel_J(idx, kb*sin(theta))) );
    }
    
    return val;  
  }
  
    
  // TODO: NEED NOT DEPEND ON GLOBAL CONSTANT
  int Quadrature_Uniform::get_N_theta( int L, double kappa, 
				       double norm_r0, double boxSize, double A, double eps )
  {
    int N_rows_max = L;
    int N_theta_max = 8*N_rows_max+1;
    //cout << "N_theta_max : " << N_theta_max << endl;
    
    // Get the transfer matrix for phi = 0, theta = [0,2PI)
    Vec3 r0(0,0,norm_r0);
    std::vector< std::vector<complex> > T = mod_transfer_fn(L, kappa, N_rows_max, 2, r0);
    
    // Copy real-space values into a single std::vector t(theta)
    std::vector<complex> t(N_theta_max);
    t[0] = T[0][0];                        // North Pole
    for( int k = 1; k < N_rows_max; ++k ) {
      t[k] = T[k][0];
      t[N_theta_max-k] = T[k][1];
    }
    
    // Get the std::absolute value of the Fourier coefficients
    ifft( &t[0], N_theta_max );
    std::vector<double> absT(N_theta_max);
    for( int k = 0; k < N_theta_max; ++k ) 
      absT[k] = std::abs(t[k]);
    
    
    // Search for N_theta
    // Want r0 = [normr0,0,0] to give same results as r0 = [0,0,norm_r0]!!
    for( int N_theta = 2*L; N_theta < N_theta_max; N_theta += 2 ) {
      
      // Get the truncation error sum_{n >= N_theta/2} + sum_{n <= -N_theta/2}
      double error_trunc = 0;
      for( int n = N_theta/2; n < N_theta_max; ++n )
	error_trunc += 2 * absT[n] * Gamma(n, std::abs(k_out)*boxSize, A);
      
      // Get the aliasing error (Only compute p=1)
      double error_alias = 2 * absT[0] * Gamma(N_theta, std::abs(k_out)*boxSize, A);
      for( int n = 1; n < N_theta/2; ++n )
	error_alias += 2 * absT[n] * Gamma(n + N_theta, std::abs(k_out)*boxSize, A);
      
      double error = 4*PI*PI/std::abs(k_out) * ( error_alias + error_trunc );
      if( error < eps ) return N_theta;
    }
    
    cerr << "Error: N_theta convergence failure" << endl;
    return N_theta_max;
  }
  
  
  int Quadrature_Uniform::get_N_phi( double& error, std::vector< std::vector<complex> >& TL, int n, 
				     int N_phiT, double maxL, double theta, int N_theta, double boxSize,
				     double A, double eps )
  {
    int N_phi;
    
    // Compute the transfer function coefficients
    ifft( &TL[n][0], N_phiT );	
    
    std::vector<double> absT(L+1);
    for( int m = 0; m <= L; ++m ) 
      absT[m] = std::abs( TL[n][m] );
    
    
    // Compute N_phi by searching for error
    for( N_phi = 4; N_phi <= ceil4( 2*maxL+1+2*LMAX ); N_phi += 4 ) {
      
      // Get the truncation error sum_{|k| >= N_phi/2}
      double error_trunc = 0;
      for( int k = N_phi/2; k <= L; ++k )
	error_trunc += 2 * absT[k] * Eta(k, std::abs(k_out)*boxSize, theta, A);
      
      // Get the aliasing error
      double error_alias = 2 * absT[0] * Eta(N_phi, std::abs(k_out)*boxSize, theta, A);
      for( int k = 1; k <= N_phi/2-1; ++k ) {
	error_alias += 2 * absT[k] * Eta(N_phi-k, std::abs(k_out)*boxSize, theta, A);
      }
      
      double error = 4*PI*PI/std::abs(k_out) * ( error_alias + error_trunc ); 
      
      if( error < eps ) break;
    }
    
    return N_phi;
  }
    
  
}

#endif
