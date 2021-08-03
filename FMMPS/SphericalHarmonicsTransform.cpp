/*
 *  SphericalHarmonicsTransform.cpp
 *  PSCAcoustic
 *
 *  Fast spherical harmonics transform.
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


#ifndef SPHERICAL_HARMONICS_TRANSFORM_CPP
#define SPHERICAL_HARMONICS_TRANSFORM_CPP

#include "Quadrature.h"
#include "SphericalHarmonicsTransform.h"

namespace FMM {
  
  
  // Constructor
  SH_Transform::SH_Transform(Quadrature* quad, Indexing* index_) : index(index_)
  {
    L = (*index)(index->size()-1, 0);
    N = quad->numRows();              // N = N_theta/2 + 1
    Q = quad->size();
    
    
    // Forward transform
    b.resize(N);
    for( int i = 0; i < N; i++ )
      b[i].resize(2*L+1);
    M.resize(N);
    
    // Double theta;
    theta.resize(N);
    Legendre.resize(N);
    int max_M = 2*L+1;
    for( int n = 0; n < N; n++ ){
      M[n] =  quad->getRow(n).numPoints();
      max_M = std::max(max_M, M[n]);
      theta[n] = quad->getRow(n).getTheta();
      Legendre[n].resize(L+1);
      
      for( int l = 0; l <= L; l++ ){
	Legendre[n][l].resize(2*l+1);
	for( int m = -l; m <= l; m++ ){
	  Legendre[n][l][m+l] = gsl_sf_legendre(l, m, cos(theta[n]));
	}
	
      }
    }
    array = new complex[max_M];
    
    // Backward transform
    c.resize((L+1)*(L+1));
    for( int i = 0; i < (L+1)*(L+1); i++ )
      c[i].resize(N);
    
    
    // ***Phi goes from 0 to 2*PI*(1-1/M[n])
    // ***Theta goes from 0 to PI
    // Phi = 2*pi*n/M[n]
    // weight = (2*pi)^2 / (M[n]*(N-1)) if not a pole
    // weight =  2*(pi)^2 / (N-1)       if pole
    weight.resize(N);
    phi.resize(N);
    for( int n = 0; n < N; n++ ){
      weight[n].assign(M[n],0.);
      phi[n].resize(M[n]);

      for( int j = 0; j < M[n]; j++ ){
	weight[n][j] = quad->getRow(n).getPoint(j).w;
	phi[n][j] = quad->getRow(n).getPoint(j).phi;
      }
    }
    
  }
  
  // Destructor
  SH_Transform::~SH_Transform()
  {
    delete array;
  }
  
  
  
  // Application of forward transform  
  std::vector<complex> SH_Transform::Forward( std::vector<complex>& vec2 )
  {
    for( int n = 0; n < N; n++ )
      for( int m = -L; m <= L; m++ )
	b[n][m+L] = 0.;
    
    assert( (int) vec2.size() == (int) index->size() );
    vec1.assign(Q, 0.);
    
    int k = 0;
    for( int n = 0; n < N; n++ ){
      for( int m = -L; m <= L; m++ ){
	for( int l = std::abs(m); l <= L; l++ ){
	  b[n][(m+2*L+1) % (2*L+1)] +=  vec2[ index->LocateIndex(l,m) ] / pow(CI, (double) (l+1)) * Legendre[n][l][m+l];
	  k++;
	}
      }
    }
    
    int idx = 0;
    for( int n = 0; n < N; n++ ){
      
      for( int i = 0; i < 2*L+1; i++ )
	array[i] = b[n][i];
      
      if( M[n] < (2*L+1) )  f_trunc(array, (2*L+1), M[n]);
      else if( M[n] > (2*L+1) ) f_extend( array, (2*L+1), M[n] );
      
      ifft(array, M[n]);
      scale(array, M[n], M[n]);
      
      
      for( int k = 0; k < M[n]; k++ ){
	vec1[idx] = array[k];
	idx++;
      }
    }
    
    return vec1;
  }
  
  
  // Application of backward transform
  std::vector<complex> SH_Transform::Backward( std::vector<complex>& vec2 ){
    assert( (int) vec2.size() == Q );
    vec1.assign(index->size(), 0.);
    
    int k = 0;
    for( int n = 0; n < N; n++ ){ 
      
      for( int j = 0; j < M[n]; j++ ){
	array[j] = weight[n][j] * vec2[k];
	k++;
      }
      
      // Perform FFT
      fft(array, M[n]);
      
      c[L][n] = array[0];
      for( int m = 1; m <= L; m++ ){
	c[L+m][n] = array[m%M[n]];
	c[L-m][n] = array[ M[n]-m%M[n] ];
      }
      
    }
    
    for( int k = 0; k < (L+1)*(L+1); k++ ){
      int l = (*index)(k,0);
      int m = (*index)(k,1);
        
      // ***Important note***
      // The real inverse requires a smoothed version of |sin()|. However, in the FMM this
      // quantity is included in the transfer function and thus should not be computed here.
      for( int n = 0; n < N; n++ )
	vec1[k] += c[L+m][n] * Legendre[n][l][m+l] * pow(CI,l+1);
    }
    
    return vec1;    
  }
  

}

#endif
