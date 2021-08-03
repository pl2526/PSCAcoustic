/*
 *  HelmholtzUtil.h
 *  PSCAcoustic
 *
 *  Routine pertaining to Helmholtz kernel.
 *
 *
 * Copyright 2014 Pierre-David Letourneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/


#ifndef HELMHOLTZUTIL_FMM_H
#define HELMHOLTZUTIL_FMM_H

#include "Vec3.h"

namespace FMM {
  
  inline complex gegenbauer_series( int L, double kappa,
				    const Vec3& r, const Vec3& r0 )
  {
    double kr0 = kappa * r0.mag();
    double kr  = kappa * r.mag();
    double rr0 = r0.dot(r) / (r0.mag() * r.mag());
    
    complex G = 0;
    int one = 1;
    for( int n = 0; n <= L; ++n ) {
      G += one * (2*n+1) * bessel_h(n,kr0) * bessel_j(n,kr) * legendre_P(n,rr0);
      one *= -1;
    }
    
    return CI*kappa * G;
  }
  
  // A function to evaluate the transfer function T_{L,r_0}(s)
  inline complex transfer_function( int L, double kappa,
				    const Vec3& r0, const Vec3& s )
  {
    assert( std::abs(s.mag()-1) < 1e-15 );
    double sdotr0 = r0.dot(s) / r0.mag();
    double knormr0 = kappa * r0.mag();
    
    complex t = 0;
    complex i_k = 1;
    for( int n = 0; n <= L; ++n ) {
      t += i_k * (2*n+1) * bessel_h(n,knormr0) * legendre_P(n,sdotr0);
      i_k *= CI;
    }
    
    return (CI*kappa)/(4*PI) * t;
  }
  
  
  // A class to quickly evaluate the transfer function T_{ell,r_0}(s) 
  // for many possible directions s
  struct Transfer_Function_Eval
  {
    int L;
    Vec3 r0hat;
    std::vector<complex> H;   // Coefficients of the series
    std::vector<double> P;
    
  Transfer_Function_Eval( complex kappa, int L_, Vec3 r0 ) 
  : L(L_), r0hat(r0/r0.mag()), H(L+1), P(L+1)
    {
      
      // TODO : My implementation does not use this factor ( (CI*kappa)/(4*PI) ) . It should be removed for efficiency
      complex i_k = (CI*kappa)/(4*PI);
      for( int n = 0; n <= L; ++n ) {
	H[n] = i_k * (2*n+1) * Amos_sf_hankel_1(n, kappa * r0.mag());  
	i_k *= CI;
      }
    }
    
    inline complex operator()(const Vec3& s)
    {
      assert( std::abs(s.mag()-1) < 1e-15 );
      double sdotr0 = r0hat.dot(s);
      gsl_sf_legendre_Pl_array(L, sdotr0, &P[0]);
      
      complex sum = 0;
      for( int n = 0; n <= L; ++n ) {
	sum += H[n] * P[n];
      }
      
      return sum;
    }
  };
  
  
  // Get the Gegenbauer truncation, L, via a number of possible methods
  inline static int get_Truncature( double kappa, 
				    double norm_r, double norm_r0, 
				    double eps, int method = 3 )
  {
    int ell = 0;
    
    switch( method ) {
    case 1: // Old Log Method
      ell = (int) (kappa*norm_r - log10(eps)*log(PI + kappa*norm_r));
      //break;
      
    case 2: // EBF Formula (Chew)
      ell = (int) (kappa*norm_r + 1.8*pow(-log10(eps), 2.0/3.0) 
		   *pow(kappa*norm_r, 1.0/3.0));
      break;
      
    case 3: // Direct method (Collino) with worst case cos = -1
      // Compute the EBF_L as an upperbound
      double knorm_r  = kappa*norm_r;
      double knorm_r0 = kappa*norm_r0;
      
      //cerr << "EBF_L = "<<get_Truncature(kappa,norm_r,norm_r0,eps,2)<<endl;
      
      // A slightly more accurate ebf_ell to use as a starting point
      //int ebf_ell = (int) (knorm_r + 1.8*pow(2-log10(eps), 2.0/3.0) 
      //	                          *pow(knorm_r, 1.0/3.0));
      
      int ebf_ell = (int) (knorm_r + 1.8*pow(-log10(eps), 2.0/3.0) 
			   *pow(knorm_r, 1.0/3.0));
      
      ell = ebf_ell;
      
      double eM = knorm_r*knorm_r0 / std::abs(norm_r0-norm_r);
      double eP = knorm_r*knorm_r0 / std::abs(norm_r0+norm_r);
      
      complex hl = bessel_h(ell, knorm_r0), hlp1 = bessel_h(ell+1, knorm_r0);
      double  jl = bessel_j(ell, knorm_r ), jlp1 = bessel_j(ell+1, knorm_r );
      
      // Decrease ell until we hit the cutoff
      double error = std::max(eM*std::abs(hlp1*jl-hl*jlp1), eP*std::abs(hlp1*jl+hl*jlp1));
      //cout << "Trunc Error: " << error << "    ell: " << ell << endl; 
      double lastError = eps;
      while( error < eps && ell > 1 ) {
	--ell;
	hlp1 = hl;  hl = bessel_h( ell, knorm_r0 );
	jlp1 = jl;  jl = bessel_j( ell, knorm_r );
	lastError = error;
	error = std::max(eM*std::abs(hlp1*jl-hl*jlp1), eP*std::abs(hlp1*jl+hl*jlp1));
	//cerr << "Trunc Error: " << error << "    ell: " << ell << endl; 
      }
      // Went one step over
      eps = lastError;
      ++ell;
      
      // Check this against the true Gegenbauer error
      // If it agrees, we're good
      // If it doesn't agree, then |r| is too small and we should
      //          do an ultradirect check (using the true Geg error)
      
      break;
    }
    
    return ell;
  }
  
  
  
  /* Computes a low-pass |sin(phi)|
   *
   * nF: number of desired frequencies of std::abs(sin)
   * nR: number of desired real space points of std::abs(sin), nR >= 2*nF+1
   */
  inline std::vector<complex> FabsSin( int nF, int nR )
  {
    assert( nR >= 2*nF-1 );
    std::vector<complex> abssin(nR,0);
    
    // Compute all the Fourier coefficients
    abssin[0] = 2.0/PI;
    for( int k = 2; k <= nF; k += 2 )
      abssin[nR - k] = abssin[k] = 2.0/(PI*(1 - k*k));
    
    fft( &abssin[0], nR );
    return abssin;
  }
  
  
  // Get the nth coefficient of |sin|
  inline double sincoef( int n ) 
  {
    if( ISODD(n) )  return 0;
    else            return 2.0/(PI*(1 - n*n));
  }
  
  
  // Takes a [-L,L] and convolves with |sin| to get [-M,M]
  inline std::vector<complex> slowConv( std::vector<complex> a, int L, int M )
  {
    ifft( &a[0], 2*L+1 );
    
    std::vector<complex> F(2*M+2,0);
    
    for( int m = 0; m <= M; ++m ) {
      for( int k = 0; k <= L; ++k ) {
	F[m] += a[k] * sincoef(m-k);
      }
      for( int k = -L; k < 0; ++k ) {
	F[m] += a[2*L+1+k] * sincoef(m-k);
      }
    }
    
    for( int m = -M; m < 0; ++m ) {
      for( int k = 0; k <= L; ++k ) {
	F[2*M+2+m] += a[k] * sincoef(m-k);
      }
      for( int k = -L; k < 0; ++k ) {
	F[2*M+2+m] += a[2*L+1+k] * sincoef(m-k);
      }
    }
    
    fft( &F[0], 2*M+2 );
    return F;
  }
  
  
  // Computes the low-pass modified transfer function matrix N_rows x N_phi
  // L:  Gegenbauer truncation
  // kappa: Wavenumber
  // N_phi: the number of quadrature points in phi
  // N_rows: the number of rows of the quadrature N_theta = 2*N_rows - 2
  // 
  // Returns: T[n][m] nth row, mth phi value of the low-pass modified trans fn
  inline std::vector< std::vector<complex> > mod_transfer_fn( int L, complex kappa,
							      int N_rows, int N_phi, 
							      const Vec3& r0 )
  {
    assert( !ISODD(N_phi) );
    
    // Construct a Transfer_Function_Eval T_{L,r_0}
    Transfer_Function_Eval Tfn( kappa, L, r0 );
    
    int N_theta = 2*N_rows - 2;           // # points used for ModTransFn
    int N_t = 2*L + 1;                    // # points needed for TransFn
    
    // Get the low-pass |sin(theta)|
    int N_sinF = (N_theta/2-1) + L;       // # freqs needed in |sin|
    int N_conv = 2*(L+N_sinF) + 1;        // size of convolution with |sin|  ??N_theta + 2L - 1??
    std::vector<complex> abssin = FabsSin(N_sinF, N_conv);
    
    double scale = 0.5 / (N_t*N_conv);    // 1/2 and interp/anterp coeffs
    
    // The N_rows x N_phi transfer function
    std::vector< std::vector<complex> > T(N_rows, std::vector<complex>(N_phi,0));
    
    // Temp storage for t_phi(theta)
    std::vector<complex> t(N_conv);
    
    // FFTs for t(theta)    (Optimized)
    fftw_complex* a = reinterpret_cast<fftw_complex*>(&t[0]);
    fftw_plan tFFT1 = fftw_plan_dft_1d(N_t,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_plan tIFFT1 = fftw_plan_dft_1d(N_conv,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_plan tFFT2 = fftw_plan_dft_1d(N_conv,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_plan tIFFT2 = fftw_plan_dft_1d(N_theta,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);
    
    
    // For each phi [0,PI)
    for( int m = 0; m < N_phi/2; ++m ) {
      double phi = m * (2*PI)/N_phi;
      double cosphi = cos( phi );
      double sinphi = sin( phi );
      
      // Sample the transfer function at 2L+1 points in theta [0,2PI)
      for( int n = 0; n < N_t; ++n ) {
	double theta = n * (2*PI)/N_t;
	double sintheta = sin( theta );
	double costheta = cos( theta );
	
	Vec3 s( sintheta*cosphi, sintheta*sinphi, costheta );
	//t[n] = scale * transfer_function(L, kappa, r0, s);
	t[n] = scale * Tfn( s );    // Optimized
      }
      
      // Interpolate the transfer function to N_conv
      fftw_execute( tFFT1 );
      f_cut( &t[0], N_t, N_conv );
      fftw_execute( tIFFT1 );
      
      // Multiply by low-pass |sin| to get the modified transfer function
      for( int k = 0; k < N_conv; ++k )
	t[k] *= abssin[k];
      
      // Anterpolate to N_theta
      fftw_execute( tFFT2 );
      f_cut( &t[0], N_conv, N_theta );
      t[N_theta/2] = 0;                    // Set the oddball freq to zero
      fftw_execute( tIFFT2 );
      
      // Explicit convolution to check the answer
      //std::vector<complex> exact = slowConv( t, L, N_theta/2-1 );
      
      // Unwrap into the N_rows x N_phi matrix
      T[0][m] = T[0][N_phi/2+m] = t[0];       // North Poles
      for( int n = 1; n < N_rows; ++n ) {
	T[n][m]         = t[n];
	T[n][N_phi/2+m] = t[N_theta-n];
      }
    }
    
    fftw_destroy_plan( tFFT1 );
    fftw_destroy_plan( tIFFT1 );
    fftw_destroy_plan( tFFT2 );
    fftw_destroy_plan( tIFFT2 );
    
    return T;
  }
  
  
}

#endif
