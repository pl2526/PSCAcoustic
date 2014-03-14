#ifndef L_CHOICE_FMM_H
#define L_CHOICE_FMM_H


#include "../Coordinates.h"
#include "TransferUtils.h"

namespace FMM {
  // TODO : Should not be here
  complex  Mie2(complex  k, complex  k_out, double radius, double rho, int l, int m = 0){
    complex  mie;
    
    complex  arg_in_C = k * radius;
    complex  arg_out_C = k_out * radius;
    complex  gamma = k/k_out * RHO_OUT/rho;
    
    // Penetrable sphere (w/ damping)
    complex  numerator = gamma * Amos_sf_bessel_jl(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C)
    - Amos_sf_bessel_jl_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C);
    complex  denominator = Amos_sf_hankel_l_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C)
    - gamma * Amos_sf_hankel_1(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C);
    
    mie = numerator / denominator;
    
    return mie;
    }
 


  
  

  
  
  int L_max_S(double r, double d, double rho, complex k, complex k_out, double eps){
    cout << "Radius : " << r << endl;
    cout << "Density : " << rho << endl;

    // TODO : What N should I use?
    // TODO : improve efficiency. statuc?
    // Indexing 
    int N = 50;
    Indexing index(N);
    
    // Representant of translation vector
    double rad, theta, phi;
    Cart2Sf(2*r+10*r, 0., 0., rad, theta, phi);

    // Compute alpha coefficients
    arma::Mat<complex> Alpha;
    Alpha = Gumerov_TransCoeff::Compute(N, rad, theta, phi, k_out, false);

    // Compute Mie coefficients
    std::vector<double> T(index.size());
    for( int i = 0; i < index.size(); i++ ){
      int l = index(i,0);
      T[i] = std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*r)) * std::abs(Mie2(k, k_out, r, rho, index(i,0), index(i,1)));
    }


    // Check error    
    // DO NOT COUNT FIRST TERM; J STARTS AT 1
    std::vector<double> err((N+1)*(N+1));
    for( int i = 0; i < (N+1)*(N+1); i++ ){
      err[i] = 0;
      for( int j = 1; j < (N+1)*(N+1); j++ ) {
	err[i] +=  T[i] * std::abs(Alpha(i,j)) / std::abs(gsl_sf_hankel_1(index(j,0), std::abs(k_out)*r)) ;
      }
    }
    
    int n = 0;
    bool flag = false;
    int L = 0;
    while( n < N ){
      
      // Check if all errors below target
      flag = true;
      for( int i = 0; i < (N+1)*(N+1); i++ ){
	if( err[i] > eps ){
	  flag = false;
	  break;
	}
      }
      
      if( flag )
	break;

      // Update L and error
      L++;
      for( int i = 0; i < (N+1)*(N+1); i++ ){
	int l = index(i,0);
	for( int m = -L; m <= L; m++ ){
	  int idx = index.LocateIndex(L,m);
	  err[i] -=  T[i] * std::abs(Alpha(i, idx)) / std::abs(gsl_sf_hankel_1(L, std::abs(k_out)*r)) ;
	}
      }
 
      n++;  
    }
    
    return L;
  };
 


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

  int L_max_T(double b, double r, double rho, complex k, complex k_out, double cte, double eps){
    

   




    // TODO : What N should I use?
    // TODO : improve efficiency. static?
    // Indexing 
    double delta = 0.8;//-6/80.*log(EPS)/log(10.);
    int N = 50;
    Indexing index(N);
    int M = (N+1)*(N+1);
    Average_H av_H(N, sqrt(3)/2.*b);
    av_H.compute();
    Average_J av_J(N, sqrt(3)/2.*b);
    av_J.compute();

    cout << "b : " << b << endl;

    double rad, theta, phi;

    Cart2Sf(0., 0., 3.5*b, rad, theta, phi);
    arma::Mat<complex> Alpha_1;
    Alpha_1 = Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), false);

    Cart2Sf(0., 0., 3.*b, rad, theta, phi);
    arma::Mat<complex> Alpha_2;
    Alpha_2 = Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), false);

    Cart2Sf(0., 0., delta*sqrt(3)/2.*b, rad, theta, phi);
    arma::Mat<complex> Beta_2;
    Beta_2 = Gumerov_TransCoeff::Compute(N, rad, theta, phi, std::abs(k_out), true);






    // Compute Mie coefficients
    std::vector<double> T(M);
    for( int i = 0; i < M; i++ ){
      int l = index(i,0);
      T[i] = std::abs(gsl_sf_hankel_1(l, std::abs(k_out)*r)) * std::abs(Mie2(k, k_out, r, rho, index(i,0), index(i,1)));
    }





    // *** Check error ***

    std::vector<double> err_1(M);
    // Initialize 1st error term
    for( int i = 0; i < M; i++ ){
	int l = index(i,0);
	int m = index(i,1);
	for( int j = 1; j < M; j++ ){
	  int l_p = index(j,0);
	  int m_p = index(j,1);

	  err_1[i] += std::abs( T[i] * Alpha_1(i,j) / gsl_sf_hankel_1(l_p, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
	}
      }



    std::vector<double> A(M);
    for( int i = 0; i < M; i++ ){
      int l = index(i,0);
      int m = index(i,1);

      A[i] = std::abs( Alpha_2(i,0) / gsl_sf_hankel_1(0, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
    }

    std::vector<double> err_2(M);
    // Initialize 2nd error term
    for( int i = 0; i < M; i++ ){
      int l = index(i,0);
      int m = index(i,1);
	
      for( int j = 1; j < M; j++ ){
	int l_p = index(j,0);
	int m_p = index(j,1);
	err_2[i] += std::abs( T[i] * Beta_2(i,j) * A[j] );
	//err_2[i] += std::abs( T[i] * Beta_2(i,j) * A[j] );
      }
    }
      




    bool flag = false;
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
       //cout << cte * err_1[i] << endl;
       //cout << cte * err_2[i] << endl;
     }

     cout << "max error : " << Max << endl;
     cout << "max error 1 : " << Max_err1 << endl;
     cout << "max error 2 : " << Max_err2 << endl;
     cout << endl;

      if( Max <= eps )
	break;


      // Update 1st error term
      for( int i = 0; i < M; i++ ){
	int l = index(i,0);
	int m = index(i,1);
	for( int m_p = -L; m_p <= L; m_p++ ){
	  int idx = index.LocateIndex(L,m_p);

	  err_1[i] -= std::abs( T[i] * Alpha_1(i,idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
	}
      }

      L++;

      // Update A
      for( int i = 0; i < M; i++ ){
	int l = index(i,0);
	int m = index(i,1);
	
	for( int m_p = -L; m_p <= L; m_p++ ){
	  int idx = index.LocateIndex(L,m_p);
	  A[i] += std::abs( Alpha_2(i,idx) / gsl_sf_hankel_1(L, std::abs(k_out)*(delta*sqrt(3.)/2.*b)));
	}
      }

      // Update 2nd error term
      for( int i = 0; i < M; i++ ){
	int l = index(i,0);
	int m = index(i,1);

	err_2[i] = 0;
	for( int j = (L+1)*(L+1); j < M; j++ ){
	  int l_p = index(j,0);
	  int m_p = index(j,1);
	  
	  err_2[i] += std::abs(T[i] * Beta_2(i,j) * A[j]);
	}
      }

      n++;  
    }
    
    return L;
  }
 





}




#endif