/*
 *  General.h
 *  PSCAcoustic
 *
 *  General library of useful methods and support packages.
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



#ifndef GENERAL_FMM_H
#define GENERAL_FMM_H

#include <cstdio>
#include <cstdlib>

#include <limits>
#include <list>
#include <algorithm>

#include <fftw3.h>

#include "../General.h"

//using namespace std;
#define list std::list



namespace FMM {
  
  // Program Constants
  const double  PI(M_PI);
  const complex CI(0,1);
  
  // Tightness of the error bound [0,1]
  //const double MLFMM_ALPHA(0.8);
  const double MLFMM_ALPHA(1);//(0.57735);

  enum Type{FORWARD, ADJOINT};
  
  //Definition Hemlholtz kernel
  struct HelmKernel
  {
    complex kappa;

    HelmKernel(){}
    ~HelmKernel(){}

    inline complex operator()( double r ) { return exp( CI * kappa * r ) / r; }
  };
  
  
#include <sys/time.h>
  
  struct StopWatch
  {
    timeval startTime, stopTime, result;
    inline void start() { gettimeofday(&startTime,NULL); }
    inline double stop() { return elapsed(); }
    inline double elapsed() {
      gettimeofday(&stopTime,NULL);
      timersub(&stopTime,&startTime,&result);
      return result.tv_sec + result.tv_usec/1000000.0;   // 10^6 uSec per Sec
    }
  };
  
  
  
  // Kill that stupid complex -op- int error
  inline complex operator*(const complex& c, const int& n) {return c*double(n);}
  inline complex operator*(const int& n, const complex& c) {return c*double(n);}
  inline complex operator/(const complex& c, const int& n) {return c/double(n);}
  inline complex operator/(const int& n, const complex& c) {return double(n)/c;}
  
#define ISODD(x) ((x) & 1)
#define ISNAN(x) ((x) != (x))
  
  // Random number in (0,1)
  inline double getRandom() {
    return drand48();
  }
  // Random number in (A,B)
  inline double getRandom( double A, double B ) {
    return A + (B-A)*getRandom();
  }
  
  inline int ceil4( double x ) {
    return 4*((int)ceil(x/4.0));
  }
  
  inline int round4( double x ) {
    return 4*((int)round(x/4.0));
  }
  
  inline int floor4( double x ) {
    return 4*((int)floor(x/4.0));
  }
  
  inline int powneg1( int n ) {
    if( ISODD(n) ) return -1;
    else           return  1;
  }
  
  
  template <class T, class Ts>
    inline T* scale( T* a, int N, Ts scale ) {
    for( int k = 0; k < N; ++k ) {
      a[k] *= scale;
    }
    return a;
  }
  
  // Overloaded complex output
  inline ostream& operator<<(ostream& os, const complex a)
  {
    //ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
    //ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
    
    int oldp = os.precision(6);
    
    os << real(a);
    if( imag(a) != 0 ) {
      if( imag(a) < 0 )
	os << " - " << -imag(a) << "*i";
      else
	os << " + " << imag(a) << "*i";
	}
    
    //os.setf(olda,ios::adjustfield);
    //os.setf(oldf,ios::floatfield);
    os.precision(oldp);
    os << "";
    
    return os;
    }
  
  
  template <class T>
    ostream& operator<<(ostream& os, const std::vector<T>& a)
    {
      ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
      ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
      
      int oldp = os.precision(8);
      
      int N = a.size();
      for( int k = 0; k < N; ++k ) {
	os << k << "\t" << a[k] << "\n";
      }
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
      os << "";
      
      return os;
    }
  
  // Shifts the zero-frequency component to the center of spectrum.
  template <class T>
    inline T* fftshift( T* a, int N ) {
    int n1 = N/2;
    T temp;
    if( ISODD(N) ) {      // N is odd
      int index = n1;
      T last = a[0];
      while( index != 0 ) {
	temp = a[index];
	a[index] = last;
	last = temp;
	index = (index + n1) % N;
      }
      a[0] = last;
    } else {          // N is even, just swaps
      for( int k = 0; k < n1; ++k ) {
	temp = a[k];
	a[k] = a[n1+k];
	a[n1+k] = temp;
      }
    }
    return a;
  }
  
  // Inverse fftshift
  template <class T>
    inline T* ifftshift( T* a, int N ) {
    int n1 = N/2;
    T temp;
    if( ISODD(N) ) {      // N is odd
      int index = ++n1;
      T last = a[0];
      while( index != 0 ) {
	temp = a[index];
	a[index] = last;
	last = temp;
	index = (index + n1) % N;
      }
      a[0] = last;
    } else {          // N is even, just swaps
      for( int k = 0; k < n1; ++k ) {
	temp = a[k];
	a[k] = a[n1+k];
	a[n1+k] = temp;
      }
    }
    return a;
  }
  
  // Truncates the sequence    N0 > N
  template <class T>
    inline T* f_trunc( T* a, int N0, int N ) {
    int M = N/2;                       // Note integer division
    int k = M + ISODD(N);
    int index = N0 - M;
    for( ; k < N; ++k, ++index ) {
      a[k] = a[index];
      a[index] = 0;
    }
    // Zero out anything that was missed
    for( ; k < N0 - M; ++k ) {
      a[k] = 0;
    }
    return a;
  }
  
  
  // Extend the sequence       N > N0
  template <class T>
    inline T* f_extend( T* a, int N0, int N ) {
    int end = N0/2 - !ISODD(N0);       // Note integer division
    int index = N0 - 1;
    int q = N - 1;
    for( ; index > end; --q, --index ) {
      a[q] = a[index];
      a[index] = 0;
    }
    // Zero out anything that was missed
    for( ; q > N0 - 1; --q ) {
      a[q] = 0;
    }
    return a;
  }
  
  template <class T>
    inline T* f_smooth( T* a, int N0, int N ) {
    if( N0 > N ) {
      return scale( f_trunc(a, N0, N), N, double(N)/N0 );
    } else if( N0 < N ) {
      return f_extend( scale(a, N0, double(N)/N0), N0, N );
    } // else they're equal and do nothing
    return a;
  }
  
  template <class T>
    inline T* f_cut( T* a, int N0, int N ) {
    if( N0 > N ) {
      return f_trunc( a, N0, N );
    } else if( N0 < N ) {
      return f_extend( a, N0, N );
    } // else they're equal and do nothing
    return a;
  }
  
  
  // Perform a (slow) forward fft
  inline complex* fft( complex* in, complex* out, int N )
  {
    fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
    fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
    fftw_plan p = fftw_plan_dft_1d( N, a, b, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute( p );
    fftw_destroy_plan( p );
    return out;
  }
  
  // Perform a (slow) in-place forward fft
  inline complex* fft( complex* in, int N ) 
  {
    return fft(in,in,N);
  }
  
  // Perform a (slow) backward fft
  inline complex* ifft( complex* in, complex* out, int N )
  {
    fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
    fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
    fftw_plan p = fftw_plan_dft_1d( N, a, b, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute( p );
    fftw_destroy_plan( p );
    scale( out, N, 1.0/N );
    return out;
  }
  
  // Perform a (slow) in-place backward fft
  inline complex* ifft( complex* in, int N )
  {
    return ifft(in,in,N);
  }
  
  inline complex* fftinterp( complex* in, int N0, int N )
  {
    if( N0 == N ) return in;
    
    fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
    fftw_plan pFFT = fftw_plan_dft_1d(N0, a, a, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pIFFT = fftw_plan_dft_1d(N, a, a, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    fftw_execute( pFFT );
    if( N0 > N )       scale( f_trunc(in, N0, N), N, 1.0/N0 );
    else if( N0 < N )  f_extend( scale(in, N0, 1.0/N0), N0, N );
    fftw_execute( pIFFT );
    
    fftw_destroy_plan( pFFT );
    fftw_destroy_plan( pIFFT );
    return in;
  }
  
  
  /* Save and Write Files */
  
  template <class T>
    inline void saveArray(T* a, int N, const char* filename) 
    {
      fstream myFile(filename, ios::out | ios::binary);
      myFile.write((char*)a, sizeof(T)*N);
      myFile.close();
    }
  
  template <class T>
    inline void saveVector( const std::vector<T> a, const char* filename )
    {
      saveArray( &a[0], a.size(), filename );
    }
  
  template <class T>
    inline T* readArray(int N, const char* filename)
    {
      T* a = new T[N];
      fstream myFile(filename, ios::in | ios::binary);
      myFile.read((char*)a, sizeof(T)*N);
      myFile.close();
      return a;
    }
  
  template <class T>
    inline std::vector<T> readVector(int N, const char* filename)
    {
      std::vector<T> a(N);
      fstream myFile(filename, ios::in | ios::binary);
      myFile.read((char*)&a[0], sizeof(T)*N);
      myFile.close();
      return a;
    }
  
  
  //////////////////
  // GSL Wrappers //
  //////////////////
  
  // Override GSL underflow error exit
  struct _GSLError_ { 
    _GSLError_() { gsl_set_error_handler(&_GSLError_::Handler); }
    static void Handler(const char* msg, const char* file, int line, int err) {
      if( err != GSL_EUNDRFLW ) {
	printf("GSLError %d in %s at %d : %s\n",err,file,line,msg);
	exit(1);
      }
    }
  };
  // Re-define GSL default error handler when loading the library
  static _GSLError_ __GSLError__;
  
  
  // Cylindrical Bessel function
  // Returns 0 in the case of underflow
  inline double bessel_J( int n, double x )
  {
    gsl_sf_result result;
    int status = gsl_sf_bessel_Jn_e(n,x,&result);
    if( status == GSL_EUNDRFLW ) return 0;
    return result.val;
  }
  
  // Spherical Bessel function j
  // Returns 0 in the case of underflow
  inline double bessel_j( int n, double x )
  {
    gsl_sf_result result;
    int status = gsl_sf_bessel_jl_e(n,x,&result);
    if( status == GSL_EUNDRFLW ) return 0;
    return result.val;
  }
  
  // Cylindral Bessel function y
  // Returns 0 in the case of underflow
  inline double bessel_Y( int n, double x )
  {
    gsl_sf_result result;
    int status = gsl_sf_bessel_Yn_e(n,x,&result);
    if( status == GSL_EUNDRFLW ) return 0;
    return result.val;
  }
  
  // Spherical Bessel function y
  // Returns 0 in the case of underflow
  inline double bessel_y( int n, double x )
  {
    gsl_sf_result result;
    int status = gsl_sf_bessel_yl_e(n,x,&result);
    if( status == GSL_EUNDRFLW ) return 0;
    return result.val;
  }
  
  // Spherical Hankel function h
  inline complex bessel_h( int n, double x )
  {
    return complex( bessel_j(n,x), bessel_y(n,x) );
  }
  
  // Legendre Polynomial
  // Returns 0 in the case of underflow
  inline double legendre_P( int n, double x )
  {
    gsl_sf_result result;
    int status = gsl_sf_legendre_Pl_e(n,x,&result);
    if( status == GSL_EUNDRFLW ) return 0;
    return result.val;
  }
  
  //Spherical harmonics function
  /*inline complex harmonic_Y( int l, int m, double theta, double phi)
  {
    gsl_sf_result result;
    int status = gsl_sf_legendre_sphPlm_e(l, std::abs(m), cos(theta), &result);
    if( status == GSL_EUNDRFLW ) return 0;
    double C = pow(-1, m);
    if( m < 0 )
      C *= pow(-1, std::abs(m));

    return C * result.val * exp(m*phi*CI);
    }*/

  
  // Gauss-Legendre quadrature
  inline void getGaussLegendreQuad(int N, std::vector<double>& x, std::vector<double>& w)
  {
    gsl_integration_glfixed_table* GLTable = 
      gsl_integration_glfixed_table_alloc(N);
    
    for( int n = 0; n < (N+1)/2; ++n ) {
      x[n + N/2] = GLTable->x[n];
      w[n + N/2] = GLTable->w[n];
    }
    for( int n = 0; n < N/2; ++n ) {
      x[n] = -x[N-1-n];
      w[n] =  w[N-1-n];
    }
    
    gsl_integration_glfixed_table_free( GLTable );
  }
  
}


#endif
