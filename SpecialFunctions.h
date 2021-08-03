/*
 *  SpecialFunctions.h
 *  PSCAcoustic
 *
 *  Special functions and more.
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


#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

#include "General.h"

// GSL library
//#include "../gsl-1.16/include/gsl/gsl_sf_bessel.h" //spherical bessel functions
#include <gsl/gsl_sf_bessel.h> //spherical bessel functions
#include "./gsl-1.16/include/gsl/gsl_sf_legendre.h" //associated Legendre polynomials (for spherical harmonics)
#include "./gsl-1.16/include/gsl/gsl_complex_math.h"
#include "./gsl-1.16/include/gsl/gsl_sf_gamma.h"
#include "./gsl-1.16/include/gsl/gsl_sf_hyperg.h"
#include "./gsl-1.16/include/gsl/gsl_sf_coupling.h"
#include "./gsl-1.16/include/gsl/gsl_integration.h"
#include "./gsl-1.16/include/gsl/gsl_complex.h"
#include "./gsl-1.16/include/gsl/gsl_sf_coupling.h"
#include "./gsl-1.16/include/gsl/gsl_math.h"
#include "./gsl-1.16/include/gsl/gsl_ieee_utils.h"
#include "./gsl-1.16/include/gsl/gsl_errno.h"
#include "./gsl-1.16/include/gsl/gsl_integration.h"

// Complex Bessel functions
//#include "../AmosBessel/libAmosBessel.h"


//GSL wrappers (Cannot handle complex  arguments)

//Hankel function of the 1st kind (gsl)
inline complex  gsl_sf_hankel_1( int l, double x){
  return gsl_sf_bessel_jl(l,x) + gsl_sf_bessel_yl(l,x) * CI;
}

inline double gsl_sf_legendre( int l, int m, double x){
  double C = 1;
  if( m < 0 )
    C *= pow(-1, std::abs(m));

  return (C * gsl_sf_legendre_sphPlm ( l, std::abs(m), x) );
}


// Y_lm (spherical harmonics)  
//\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) exp(I*m*phi) )
inline complex  gsl_sf_harmonic( int l, int m, double theta, double phi){
  if( l >= std::abs(m) ){
    return gsl_sf_legendre(l, m, cos(theta)) * exp(m*phi*CI);
  } else {
    return 0.;
  }
}

inline complex  gsl_sf_clebsch(int l1, int l2, int l3, int m1, int m2, int m3){
  return pow(-1., (double) l1-l2+m3) * sqrt(2.*l3+1.) * gsl_sf_coupling_3j(2*l1, 2*l2, 2*l3, 2*m1, 2*m2, -2*m3);
}

inline double gsl_sf_gaunt(int n, int m, int nu, int mu, int q){
  return pow(-1., (double) (m+mu) ) * sqrt((2*n+1)*(2*nu+1)*(2*q+1)/(4*PI)) * 
    gsl_sf_coupling_3j(2.*n,2.*nu,2.*q,0,0,0) * gsl_sf_coupling_3j(2.*n,2.*nu,2.*q,2.*m,2.*mu,-2.*(m+mu));
}


// TODO: FIX AMOS
//Amos wrappers for Bessel functions (these can handle complex  arguments)
inline complex  Amos_sf_bessel_jl(int l, complex  z){
  	complex  value;
	//AmosBessel('j', z, (double) l, 1, 0, &value);
	value = gsl_sf_bessel_jl(l,std::real(z));
	return value;
}

inline complex  Amos_sf_bessel_yl(int l, complex  z){
  	complex  value;
	//AmosBessel('y', z, (double) l, 1, 0, &value);
	value = gsl_sf_bessel_yl(l,std::real(z));
	return value;
}


inline complex  Amos_sf_hankel_1(int l, complex  z){
  //complex  valueR;
	//AmosBessel('j', z, (double) l, 1, 0, &valueR);
  	//complex  valueI;
	//AmosBessel('y', z, (double) l, 1, 0, &valueI);
	//return (valueR + valueI*CI);

	complex value = gsl_sf_hankel_1(l, std::real(z));
	return value;
}


// Derivatives

inline double gsl_sf_bessel_jl_prime(int l, double z){
  return l/z * gsl_sf_bessel_jl(l,z) - gsl_sf_bessel_jl(l+1,z);
}

inline complex  gsl_sf_hankel_l_prime(int l, double z){
  return l/z * gsl_sf_hankel_1( l, z) - gsl_sf_hankel_1(l+1, z);
}

inline complex  Amos_sf_bessel_jl_prime(double l, complex  z){
  return l/z * Amos_sf_bessel_jl(l,z) - Amos_sf_bessel_jl(l+1,z);
}

inline complex  Amos_sf_hankel_l_prime(double l, complex  z){
  return l/z * Amos_sf_hankel_1( l, z) - Amos_sf_hankel_1(l+1, z);
}




// Chebyshev polynomial evaluation
inline double Cheb(int n, double x){
  return cos(n*acos(x));
}







// Gauss-Legendre quadrature
inline void getGaussLegendreQuad(int N, std::vector<double>& x, std::vector<double>& w)
{
  x.resize(N);
  w.resize(N);

  gsl_integration_glfixed_table* GLTable = gsl_integration_glfixed_table_alloc(N);
  
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

// Gauss-Chebyshev quadrature
inline void getGaussChebyshevQuad(int N, std::vector<double>& x, std::vector<double>& w)
{
  x.resize(N);
  w.resize(N);

  for( int n = 0; n < N; n++ ) {
    x[n] = cos( (2.*n+1.)/(2.*N)*PI  );
    w[n] = 2./N;
  }
}




//---------General Functions--------//

#define ISODD(x) ((x) & 1)
#define ISNAN(x) ((x) != (x))




//TODO : still needed?
inline int choose( int n, int r ){
  
  if( (n < r) || (r < 0) || (n < 0) ){
    return 0;
  } else {
    return gsl_sf_gamma(n+1) / (gsl_sf_gamma(r+1) * gsl_sf_gamma(n-r+1));
  }
  
}






// Implementation os Sparse Matrix class
// TODO:  this takes a lot of memory?
class SpMatrix
{
  map<std::pair<int,int>, complex > M;

 public:

  // Empty constructor
  SpMatrix(){}

  ~SpMatrix(){}

  void add(int i, int j, complex  val){
    assert( i >= 0 );
    assert( j >= 0 );

    std::pair<int,int> coord = std::make_pair(i,j);
    if( M.count(coord) ){
      M.insert(std::make_pair(coord, val));
    } else {
      M[coord] = val;
    }

  }


  
  complex  operator () (int i, int j){
    std::pair<int,int> coord = std::make_pair(i,j);
    complex  val = 0.;

    ( M.count(coord) ) ? val = M.find(coord)->second : val = 0.;
    
    return val;
  }

  inline int size() { return M.size(); }

};




#endif

