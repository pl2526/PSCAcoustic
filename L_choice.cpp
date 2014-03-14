/*
 *  L_choice.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 12/3/12.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Routines used to established the number of coefficients necessary for each scatterers
 *  and at each level in order to achieve a certain accuracy
 */
#ifndef L_CHOICE_CPP
#define L_CHOICE_CPP

#include "L_choice.h"

// ***Auxilliaries***

double ABS_J (double x, void * params) {
  int l = *(int *) params;
  return  gsl_sf_bessel_jl(l, x);
}

struct Average_J{
  
  std::vector<double> int_J;
  int N;
  double R;    
  
  Average_J( int N_, double R_):N(N_), R(R_){}
  
  void compute(){
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(5000);
    
    int_J.resize(N+1);
    
    double val;
    double err;
    
    double abserr = 1e-15;
    double relerr = 1e-12;
    
    gsl_function F;
    F.function = &ABS_J;
    
    for( int l = 0; l <= N; l++ ){
      F.params = &l;
      
      gsl_integration_qag(&F, 0, R, abserr, relerr, 5000, 1, workspace, &val, &err);
      int_J[l] = val;
      
      //cout << val << endl;
    }
    
    gsl_integration_workspace_free(workspace);
  }
  
  inline double get(int l, int l_p){
    //double val = 0;
    //for( int l_pp = std::abs(l-l_p); l_pp <= (l+l_p); l_pp++ )
    //val += int_J[l_pp];
    
    return  int_J[l_p];
  }
  
};



double RECIPROC_H (double x, void * params) {
  int l = *(int *) params;
  return  4*PI*x*x / std::abs(gsl_sf_hankel_1(l, x));
}

struct Average_H{
  
  std::vector<double> int_H;
  int N;
  double R;    
  
  Average_H( int N_, double R_):N(N_), R(R_){}
  
  void compute(){
    
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(5000);
    
    int_H.resize(N+1);
    
    double val;
    double err;
    
    double abserr = 1e-15;
    double relerr = 1e-12;
    
    gsl_function F;
    F.function = &RECIPROC_H;
    
    for( int l = 0; l <= N; l++ ){
      F.params = &l;
      
      gsl_integration_qag(&F, 0, R, abserr, relerr, 5000, 1, workspace, &val, &err);
      int_H[l] = std::abs(gsl_sf_hankel_1(l, R)) * val;
      
      //cout << val << endl;
    }
    
    gsl_integration_workspace_free(workspace);
  }
  
  inline double get(int l){
    return int_H[l];
  }
  
};



double av_beta(int l, int l_p, int m, int m_p, double r){
  double val = 0;
  if( l == l_p && m == m_p )
    val = 4*PI * ( sin(r) - cos(r) );
  
  return val;
}

double av_alpha(int l, int l_p, int m, int m_p, double r, double R){
  double val = 0;
  if( l == l_p && m == m_p )
    val = 4*PI * std::abs( exp(r*CI)*(r+CI) - exp(R*CI)*(R+CI) );
  
  return val;
}




// ***Main functions***

/*
  int L_max_scat(double r, double d, double C, double A, int N = 50){
  
  std::vector<double> W(N);
  for( int l = 0; l < N; l++ )
  W[l] = 1. / std::abs(Amos_sf_hankel_1(l, K_L_OUT*r));


  int L = 0;
  bool precision;
  double bound;
  double val;
  double Delta = std::max(std::real(K_L_OUT),std::real(K_T_OUT))*(2.*r+d);

  Pvec p(0.,0.,0.,Pvec::Cartesian);
  Scatterer S(RADIUS, K_L, K_T, K_L_OUT, K_T_OUT, RHO, p);
  vector< std::vector<complex> > TM S.getTMatrix();

  while( L < N ){
    precision = true;

    for( int l_s = 0; l_s < N; l_s++ ){
	bound = 0;

	val = 1;
	int l = L+1;


	double T =  4*PI*C*std::abs(Mie(K_L, K_L_OUT, r, RHO, l_s));
	while( val > 1e-20 || l < N ){


	  val = 0;
	  //for( int m = -l; m <= l; m++ ){
	    if( (l_s <= L && l > L) || (l_s > L) ){
	      double h = T * W[l] * H_bound2(l, l_s, Delta);
	      bound += h;
	      val = std::max(val, std::abs(h));
	    }

	  l++;
	}
	  

	if( A*bound > EPS ){
	  precision = false;
	  break;
	}

    }
    
    if( precision )
      break;
    else
      L++;

    cout << L << endl;
  }
  
  return L;
}
*/



int L_max_level(double b, double cte, complex k_out, int N ){
    

    // TODO : What N should I use?
    // TODO : improve efficiency. static?
    // Indexing 
    double delta = 0.8;//-6/80.*log(EPS)/log(10.);
    Indexing index(N);

    // TODO: Why squared?
    int M = (N+1)*(N+1);

    Average_H av_H(N, sqrt(3)/2.*b);
    av_H.compute();
    Average_J av_J(N, sqrt(3)/2.*b);
    av_J.compute();


    double rad, theta, phi;

    Cart2Sf(0., 0., 3.5*b, rad, theta, phi);
    arma::Mat<complex> Alpha_1;
    Alpha_1 = FMM::Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), false);

    Cart2Sf(0., 0., 3.*b, rad, theta, phi);
    arma::Mat<complex> Alpha_2;
    Alpha_2 = FMM::Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), false);

    Cart2Sf(0., 0., delta*sqrt(3)/2.*b, rad, theta, phi);
    arma::Mat<complex> Beta_2;
    Beta_2 = FMM::Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), true);






    // Compute bound on T-matrix coefficients
    std::vector<complex> T(M);
    Scatterer::TMconstruct(T, K, K_OUT, RADIUS, RHO);
    for( int i = 0; i < M; i++ ){
      //int l = index(i,0);
      T[i] = std::abs(T[i]);
      //T[i] = std::max(T[i], std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*RADIUS) * TM[i][j]));
      
      //T[i] = std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*r)) * std::abs(Mie2(k, k_out, r, rho, index(i,0), index(i,1)));
    }
    



    // *** Check error ***

    std::vector<double> err_1(M);
    // Initialize 1st error term
    for( int i = 0; i < M; i++ ){
      for( int j = 1; j < M; j++ ){
	int l_p = index(j,0);
	err_1[i] += std::abs( T[i] * Alpha_1(i,j) / gsl_sf_hankel_1(l_p, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
      }
    }

    

    std::vector<double> A(M);
    for( int i = 0; i < M; i++ )
      A[i] = std::abs( Alpha_2(i,0) / gsl_sf_hankel_1(0, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));

    std::vector<double> err_2(M);
    // Initialize 2nd error term
    for( int i = 0; i < M; i++ )
      for( int j = 1; j < M; j++ )
	err_2[i] += std::abs( T[i] * Beta_2(i,j) * A[j] );




    int L = 0;
    int n = 0;
    while( n < N ){

     double Max = 0;
     double Max_err1 = 0;
     double Max_err2 = 0;
     for( int i = 0; i < M; i++ ){
       Max = std::max(Max, cte * ( err_1[i] + err_2[i] ));
       //Max = std::max(Max, cte * ( err_1[i] + err_2[i]) );
       Max_err1 = std::max(Max_err1, cte * err_1[i]);
       Max_err2 = std::max(Max_err2, cte * err_2[i]);
     }

     cout << "max error : " << Max << endl;
     //cout << "max error 1 : " << Max_err1 << endl;
     //cout << "max error 2 : " << Max_err2 << endl;
     //cout << endl;

     // Check if stopping criterion is satisfied
      if( Max <= EPS )
	break;


      // Update 1st error term
      for( int i = 0; i < M; i++ ){
	for( int m_p = -L; m_p <= L; m_p++ ){
	  int idx = index.LocateIndex(L,m_p);
	  err_1[i] -= std::abs( T[i] * Alpha_1(i,idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
	}
      }

      L++;

      // Update A
      for( int i = 0; i < M; i++ ){
	for( int m_p = -L; m_p <= L; m_p++ ){
	  int idx = index.LocateIndex(L,m_p);
	  A[i] += std::abs( Alpha_2(i,idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
	}
      }

      // Update 2nd error term
      for( int i = 0; i < M; i++ ){
	err_2[i] = 0;
	for( int j = (L+1)*(L+1); j < M; j++ )
	  err_2[i] += std::abs(T[i] * Beta_2(i,j) * A[j]);
      }

      n++;  
    }
    
    return L;
  }




#endif
