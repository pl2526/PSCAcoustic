#ifndef RECTEST_CPP
#define RECTEST_CPP

/*
 *  MAIN.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Overstructure.
 */

#include <vector>
#include <complex>
#include <iostream>

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "PSCEnv.h"
#include "ElasticCoordChange.h"
#include "Solver.h"
#include "Finalize.h"

Indexing temp_index(LMAX_VSWF);

// p: point at which to evaluate expansion
// c: center of the expansion
std::vector<complex> Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a_0, std::vector< std::vector<complex> >& a_1, 
			      std::vector< std::vector<complex> >& a_2, ECC::type Type){
  std::vector<complex> v(3); v[0] = 0.; v[1] = 0.; v[2] = 0.;
  std::vector< std::vector<complex> > curl;
  std::vector< std::vector<complex> > curl_curl;
  std::vector< std::vector<complex> > grad;

  int N = a_0.size();
  for( int n = 0; n < N; n++ ){
    Pvec q = p-ScL[n].location;
    curl = ECC::curl_VSWF(temp_index, q, K_L_OUT, Type);
    curl_curl = ECC::curl_curl_VSWF(temp_index, q, K_L_OUT, Type);
    grad = ECC::gradient_VSWF(temp_index, q, K_L_OUT, Type);

    for( int i = 0; i < (LMAX_VSWF+1)*(LMAX_VSWF+1); i++ ){
      v[0] += curl[i][0]*a_0[n][i] + curl_curl[i][0]*a_1[n][i] + grad[i][0]*a_2[n][i];
      v[1] += curl[i][1]*a_0[n][i] + curl_curl[i][1]*a_1[n][i] + grad[i][1]*a_2[n][i];
      v[2] += curl[i][2]*a_0[n][i] + curl_curl[i][2]*a_1[n][i] + grad[i][2]*a_2[n][i];
    }
  }

  return v;
}


int main(int argc,char **args)
{

  // Display relevant information concerning the problem
  cout << endl << endl;
  cout << "RTOL : " << RTOL << endl;
  cout << "MAXITS : " << MAXITS << endl;
  cout << "OMEGA : " << OMEGA << endl;
  cout << "LAMBDA_OUT : " << LAMBDA_OUT << endl;
  cout << "MU_OUT : " << MU_OUT << endl;
  cout << "C_L_OUT : " << C_L_OUT << endl;
  cout << "C_T_OUT : " << C_T_OUT << endl;
  cout << "K_L_OUT : " << K_L_OUT << endl;
  cout << "K_T_OUT : " << K_T_OUT << endl;
  cout << "NScat : " << NScat << endl;
  cout << "RADIUS : " << RADIUS << endl;
  cout << "nLevels : " << nLevels << endl;
  cout << "EPS : " << EPS << endl;

  cout << endl;


  // PSC environment for fast algorithm
  PSC_Env PSC_env(nLevels, EPS); 

  // Linear solver
  Solver solver; 


  // TODO: Should not have different instances of Index
  Indexing Index(LMAX_VSWF);


  // Construct scatterers
  cout << "***Building scatterers..." << endl;
  std::vector<Scatterer> ScL_A = ScatInit( NScat, K_L, K_T, K_L_OUT, K_T_OUT );
  std::vector<Scatterer> ScL_B = ScatInit( NScat, K_L, K_T, K_L_OUT, K_T_OUT );
  cout << "***Building scatterers: done" << endl << endl;


  
  // FIRST FIELD

  // Construct PSC environment
  cout << "***Constructing environment..." << endl;
  PSC_env.Construct(K_L_OUT, K_T_OUT, ScL_A);
  cout << "***Constructing environment: done" << endl << endl;
 

  cout << "***Computing solution..." << endl;
  
  cout << "   ***Initializing source (plane wave)..." << endl;
  // Initialization of source (plane waves)
  std::string src_type("PW");
  Pvec direction(1./sqrt(2.), 1./sqrt(2.), 0.,Pvec::Cartesian);                // Direction of propagation
  double long_comp = 1.;
  Pvec perp_comp(-1., 1., 0.5, Pvec::Cartesian);          // Shear component perpendicular to direction of propagation 
  IncWave* IW_A = new PlaneWave(solver.getECC(), K_L_OUT, K_T_OUT, direction, long_comp, perp_comp);   
  cout << "   ***Initializing source: done" << endl << endl;
  
  
  // Initialization of right-hand side
  cout << "   ***Initializing right-hand side..." << endl;
  std::vector< std::vector<complex > > RHS_0_A(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for curl VSWF
  std::vector< std::vector<complex > > RHS_1_A(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for curl-curl VSWF
  std::vector< std::vector<complex > > RHS_2_A(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for gradient VSWF

  solver.RHSInit(IW_A, ScL_A, RHS_0_A, RHS_1_A, RHS_2_A);  
  cout << "   ***Initializing right-hand side: done" << endl << endl;
  
  
  // Vectors storing solutions for each kind of wave
  // Convention: 0 -> curl VSWF; 1 -> curl_curl VSWF; 2 -> gradient VSWF; 
  cout << "   ***Solving linear system..." << endl;
  std::vector< std::vector<complex> > u_0_A(NScat, std::vector<complex>(Index.size()));
  std::vector< std::vector<complex> > u_1_A(NScat, std::vector<complex>(Index.size()));
  std::vector< std::vector<complex> > u_2_A(NScat, std::vector<complex>(Index.size()));

  // TODO: send solve to a file
  double res_A = 1e10, rel_res_A = 1e10, cond;
  solver.Solve(PSC_env, ScL_A, RHS_0_A, RHS_1_A, RHS_2_A, u_0_A, u_1_A, u_2_A, res_A, rel_res_A, cond);
  cout << "   ***Solving linear system: done" << endl << endl;
  
  
  // Write information about problem to file
  /* if( rank == 0){
     cout << "   ***Writing to file..." << endl;
     
     Write_Info(Proc_Idx, src_type, direction, long_comp, perp_comp, res, rel_res, cond);
     Write_Location(Proc_Idx, ScL);
     Write_Solution( Proc_Idx, Index, u_0, u_1, u_2); // Solution OUTSIDE scatterers
     
     cout << "   ***Writing to file: done" << endl << endl;
     } */
  

  // SECOND FIELD

 // Construct PSC environment
  cout << "***Constructing environment..." << endl;
  PSC_env.Destruct();
  PSC_env.Construct(K_L_OUT, K_T_OUT, ScL_B);
  cout << "***Constructing environment: done" << endl << endl;

  cout << "***Computing solution..." << endl;
  
  cout << "   ***Initializing source (plane wave)..." << endl;
  // Initialization of source (plane waves)
  Pvec direction_B(-1./sqrt(2.), 0., 1./sqrt(2.),Pvec::Cartesian);                // Direction of propagation
  double long_comp_B = 0.25;
  Pvec perp_comp_B(-2., 0.25, 2, Pvec::Cartesian);          // Shear component perpendicular to direction of propagation 
  IncWave* IW_B = new PlaneWave(solver.getECC(), K_L_OUT, K_T_OUT, direction_B, long_comp_B, perp_comp_B);   
  cout << "   ***Initializing source: done" << endl << endl;
  
  
  // Initialization of right-hand side
  cout << "   ***Initializing right-hand side..." << endl;
  std::vector< std::vector<complex > > RHS_0_B(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for curl VSWF
  std::vector< std::vector<complex > > RHS_1_B(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for curl-curl VSWF
  std::vector< std::vector<complex > > RHS_2_B(NScat, std::vector<complex >(Index.size()));    // Multipole coefficients for gradient VSWF
  
  solver.RHSInit(IW_B, ScL_B, RHS_0_B, RHS_1_B, RHS_2_B);  
  cout << "   ***Initializing right-hand side: done" << endl << endl;
  
  
  // Vectors storing solutions for each kind of wave
  // Convention: 0 -> curl VSWF; 1 -> curl_curl VSWF; 2 -> gradient VSWF; 
  cout << "   ***Solving linear system..." << endl;
  std::vector< std::vector<complex> > u_0_B(NScat, std::vector<complex>(Index.size()));
  std::vector< std::vector<complex> > u_1_B(NScat, std::vector<complex>(Index.size()));
  std::vector< std::vector<complex> > u_2_B(NScat, std::vector<complex>(Index.size()));
  
  double res_B = 1e10, rel_res_B = 1e10, cond_B;
  solver.Solve(PSC_env, ScL_B, RHS_0_B, RHS_1_B, RHS_2_B, u_0_B, u_1_B, u_2_B, res_B, rel_res_B, cond_B);
  cout << "   ***Solving linear system: done" << endl << endl;
  





  // POST-PROCESSING
  
  // Create a grid onto which to compute field
  double R = 1;
  int M = 1e2;
  double delta = pow(2.*R/M,3.);
  std::vector<double> X(2*M+1);
  for( int i = -M; i < M; i++ )
       X[i] = R * i/M;

  // Compute volume integral
  complex int_V = 0.;
  complex int_S = 0.;

  std::vector<complex> v_pw_A;
  std::vector<complex> v_pw_B;
  std::vector<complex> v_sc_A;
  std::vector<complex> v_sc_B;

  complex dot1;
  complex dot2;

  for( int i = 0; i < M; i++ ){
    cout << "i : " << i << endl;
    for( int j = 0; j < M; j++ ){
      for( int k = 0; k < 1; k++ ){
	if( sqrt(X[i]*X[i] + X[j]*X[j] + X[k]*X[k]) < R ){
	  Pvec p(X[i],X[j],X[k],Pvec::Cartesian);
	  
	  // Evaluate incoming field at point
	  v_pw_A = IW_A->Evaluate(p);
	  v_pw_B = IW_A->Evaluate(p);
	  
	  // Evaluate scattered field at point
	  v_sc_A = Evaluate(p, ScL_A, u_0_A, u_1_A, u_2_A, ECC::Hankel);
	  v_sc_B = Evaluate(p, ScL_B, u_0_B, u_1_B, u_2_B, ECC::Hankel);
	  
	  //cout << v_pw_A[0] << " : " <<v_pw_A[1] << " : " <<v_pw_A[2] << endl;
	  //cout << v_pw_B[0] << " : " <<v_pw_B[1] << " : " <<v_pw_B[2] << endl << endl;
	  
	  //cout << v_sc_A[0] << " : " <<v_sc_A[1] << " : " <<v_sc_A[2] << endl;
	  //cout << v_sc_B[0] << " : " <<v_sc_B[1] << " : " <<v_sc_B[2] << endl << endl;
	  
	  // Compute dot products
	  dot1 = 2*RHO_OUT*( v_sc_A[0]*std::conj(v_pw_B[0]) + v_sc_A[1]*std::conj(v_pw_B[1]) + v_sc_A[2]*std::conj(v_pw_B[2]) );
	  dot2 = 2*RHO_OUT*( v_sc_B[0]*std::conj(v_pw_A[0]) + v_sc_B[1]*std::conj(v_pw_A[1]) + v_sc_B[2]*std::conj(v_pw_A[2]) );
	  
	  int_V = delta * (dot1 + dot2);
	}
      }
    }
  }

  cout << "int_V : " << int_V << endl;
  






  
  delete IW_A;
  delete IW_B;

  // Clear memory
  PSC_env.Destruct();


  return 0;
}






#endif
