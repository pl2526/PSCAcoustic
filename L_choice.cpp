/*
 *  L_choice.cpp
 *  PSCAcoustic
 *
 *  Routines used to established the number of coefficients necessary for each scatterers
 *  and at each level in order to achieve a certain accuracy
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

#ifndef L_CHOICE_CPP
#define L_CHOICE_CPP

#include "L_choice.h"


int L_max_scat(double boxSize, double d_min, double eps, complex k, complex k_out, int N){

  double beta = 0.25;
  double delta = 2.5;

  N = 50;
  std::vector<double> err((N+1)*(N+1));
  FMM::Z_transfer z_transfer;
  SpMatrix Alpha_Z;
  z_transfer.Compute(Alpha_Z, N, delta*boxSize, K_OUT, false);

  std::vector<double> T(N);
  for( int l = 0; l < N; l++ )
    T[l] = std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*(RADIUS + beta*boxSize) )) 
      * std::abs(Scatterer::Mie(k, k_out, RADIUS, RHO, 0, l, 0));

  int i = 0;
  for( int l = 0; l < N; l++ ){
    for( int m = -l; m <= l; m++ ){
      err[i] = 0;
      
	 int j = 0;
	 for( int l_p = 0; l_p < N; l_p++ ){
	   for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	     if( m == m_p )
	       err[i] += std::abs( T[l] * Alpha_Z(i,j) / Amos_sf_hankel_1(l_p, std::abs(k_out)*(RADIUS + beta*boxSize) ) ) ;
	     j++;
	   }
	 }
	 
	 i++;
    }
  }
  

    int L = 0;
    while( L < N ){
      
      double Max = 0;
      for( int i = 0; i < (N+1)*(N+1); i++ )
	Max = std::max(Max, err[i]);
      cout << "Max : " << Max << endl;
      
      if( Max <= eps )
	break;
      
      
      int i = 0;
      for( int l = 0; l < N; l++ ){
	for( int m = -l; m <= l; m++ ){
	  
	  for( int m_p = -L; m_p <= L; m_p++ ){
	    if( m == m_p ){
	      int idx = L*L + (L + m_p);
	      err[i] -= std::abs( T[l] * Alpha_Z(i, idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(RADIUS + beta*boxSize)) );
	      }
	  }

	  i++;
	}
      }
      
      L++;
    }
    
    return L;
}



// TODO: Pass Mie flag
int L_max_level(double boxSize, double eps, double radius, complex k, complex k_out, int N ){

  cout << "boxSize : " << boxSize << endl;
 cout << "eps : " << eps << endl;
 cout << "k: " << k << endl;
 cout << "k_out : " << k_out << endl;
 cout << "N : " << N << endl;
    
 double beta_in = 1.1;
 double beta_out = 1.1;//0.05/radius;
 double beta = 0.05;
    double delta = 0.75;
    N =  std::max(N + 0., 10*std::ceil(std::abs(k_out)/2/PI*boxSize));

    FMM::Z_transfer z_transfer;
    SpMatrix Alpha_1_Z;
    z_transfer.Compute(Alpha_1_Z, N, 3.5*boxSize, K_OUT, false);
    cout << "foo1" << endl;

    SpMatrix Alpha_2_Z;
    z_transfer.Compute(Alpha_2_Z, N, 3.*boxSize, K_OUT, false);
    cout << "foo2" << endl;
    //Alpha_1_Z = Alpha_2_Z;

    SpMatrix Beta_2_Z;
    z_transfer.Compute(Beta_2_Z, N, delta*boxSize, K_OUT, true);
    cout << "foo3" << endl;



    // Pre-Compute Mie coefficients
    std::vector<double> T(N);
    for( int l = 0; l < N; l++ )
      T[l] = std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*(radius+beta*boxSize) )) * std::abs(Scatterer::Mie(k, k_out, radius, RHO, 0, l, 0));
    //T[l] = std::abs(gsl_sf_hankel_1(l, beta_out*std::abs(k_out)*radius)) * std::abs(Scatterer::Mie(k, k_out, radius, RHO, 0, l, 0));




    // *** Check error ***
    int M = (N+1)*(N+1);
    std::vector<double> err_1(M);
    std::vector<double> A(M);
    std::vector<double> err_2(M);

    int i = 0;
    for( int l = 0; l < N; l++ ){
      for( int m = -l; m <= l; m++ ){

	int j = 1;
	for( int l_p = 1; l_p < N; l_p++ ){
	  for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	    if( m == m_p )
	      err_1[i] += std::abs( T[l] * Alpha_1_Z(i, j) / gsl_sf_hankel_1(l_p, std::abs(k_out)*(radius+beta*boxSize) ) );
	    //err_1[i] += std::abs( T[l] * Alpha_1_Z(i, j) / gsl_sf_hankel_1(l_p, beta_in*std::abs(k_out)*radius) );
	    j++;
	  }
	}

	i++;
      }
    }



     i = 0;
      for( int l = 0; l < N; l++ ){
	for( int m = -l; m <= l; m++ ){
	  A[i] = 0;

	  int j = 0;
	  for( int l_p = 0; l_p < N; l_p++ ){
	    for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	      if( m == m_p )
		A[i] += std::abs( Alpha_2_Z(i,j) / gsl_sf_hankel_1(l_p, delta*std::abs(k_out)*boxSize));
	      j++;
	    }
	  }

	  i++;
	}
      }

     i = 0;
     for( int l = 0; l < N; l++ ){
       for( int m = -l; m <= l; m++ ){
	 err_2[i] = 0;

	 //int j = L*L+2*L+1;
	 int j = 0;
	 for( int l_p = 0; l_p < N; l_p++ ){
	   for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	     if( m == m_p )
	       err_2[i] += std::abs(T[l] * Beta_2_Z(i,j) * A[j]);
	     j++;
	   }
	 }
	 
	 i++;
       }
     }
     



    bool flag = false;
    int L = 0;
    while( L < N ){


     double Max = 0;
     double Max_err1 = 0;
     double Max_err2 = 0;
     for( int i = 0; i < M; i++ ){
	Max = std::max(Max, err_1[i] + err_2[i]);
	Max_err1 = std::max(Max_err1, err_1[i]);
	Max_err2 = std::max(Max_err2, err_2[i]);
     }

 
     if( Max <= eps && L >= 1 )
	break;

      int i = 0;
      for( int l = 0; l < N; l++ ){
	for( int m = -l; m <= l; m++ ){

	  for( int m_p = -L; m_p <= L; m_p++ ){
	    if( m == m_p ){
	      int idx = L*L + (L + m_p);
	      err_1[i] -= std::abs( T[l] * Alpha_1_Z(i, idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(radius+beta*boxSize) ) );
	      //err_1[i] -= std::abs( T[l] * Alpha_1_Z(i, idx) / gsl_sf_hankel_1(L, beta_in*std::abs(k_out)*radius) );
	      }
	  }

	  i++;
	}
      }




      i = 0;
      for( int l = 0; l < N; l++ ){
	for( int m = -l; m <= l; m++ ){

	    for( int m_p = -L; m_p <= L; m_p++ ){
	      if( m == m_p){
		int j = L*L + L + m_p;
		err_2[i] -= std::abs(T[l] * Beta_2_Z(i,j) * A[j]);
	      }
	    }

	  i++;
	}
      }


      L++;
    }
    
    return L;
}




#endif
