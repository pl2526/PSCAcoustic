/*
 *  SphericalHarmonicsTransform.cpp
 *  PSCAcoustic
 *
 *  Fast spherical harmonics transform.
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
