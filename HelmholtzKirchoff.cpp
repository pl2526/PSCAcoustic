#ifndef HK_CPP
#define HK_CPP

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
#include "Solver.h"
#include "Finalize.h"

Indexing temp_index(LMAX);



int main(int argc,char **args)
{
  int Proc_Idx = 0;

  // Properties of sampling manifold
  double R = 1e9;     // Radius of spherical manifold
  int size_t = 500;   // Number of smaple points in theta
  int size_p = 500;   // Number of smaple points in phi

  // Properties of source
  std::string src_type("PtW");
  //double amplitude = 1;
  Pvec src_location(25., 30., 40.,Pvec::Cartesian);
  complex src_amplitude = 1.;


  // Properties of images
  double L = 30;    // Size of window [-L,L]x[-L,L]
  int M = 150;         // Number of samples per dimension (2*M+1)

  // Display relevant information concerning the problem
  cout << endl << endl;
  cout << "LMAX : " << LMAX << endl;
  cout << "RTOL : " << RTOL << endl;
  cout << "MAXITS : " << MAXITS << endl;
  cout << "OMEGA : " << OMEGA << endl;
  cout << "C_OUT : " << C_OUT << endl;
  cout << "C : " << C << endl;
  cout << "K_OUT : " << K_OUT << endl;
  cout << "K : " << K << endl;
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
  Indexing Index(LMAX);

  // Construct scatterers
  cout << "***Building scatterers..." ;
  std::vector<Scatterer> ScL;
  Pvec center(0.,0.,0. ,Pvec::Cartesian);                                          // Center of cluster
  std::vector<Pvec> scat_loc = RandSpherical(20., RADIUS, D_MIN, center, NScat);    // Random distribution
 assert( scat_loc.size() == NScat );
 
 // Construct scatterers
 for( int n = 0; n < NScat; n++ ){
   Scatterer scatterer(RADIUS, K, K_OUT, RHO, scat_loc[n]);
   ScL.push_back(scatterer);
 }
 cout << "* Done" << endl << endl;

  // Construct PSC environment
  cout << "***Constructing environment..." << endl;
  PSC_env.Construct(K_OUT, ScL);
  cout << "***Constructing environment: done" << endl << endl;
 
  



  // *** Solve forward problem


  // ----For all that follows:
  // P is used to denote computations associated with the pressure Green's function, i.e., G^P (Ammari et al.) 
  // S is used to denote computations associated with the shear Green's function, i.e., G^S (Ammari et al.) 


  // Initialization of source
  // R.h.s.
  SphericalWave* IW = new SphericalWave(K_OUT, src_amplitude, src_location);
  std::vector< std::vector<complex > > RHS(NScat, std::vector<complex >(Index.size())); 
  solver.RHSInit(IW, ScL, RHS);  


  
  // Compute G
  cout << "   ***Solving linear system..." << endl;
  std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));

  // Solve linear system
  cout << " Computing solution forward problem..." << endl;
  double res = 1e10, rel_res = 1e10, cond;
  solver.Solve(PSC_env, ScL, RHS, u, res, rel_res, cond);
  cout << " Computing solution P forward problem: done" << endl;
  
  // Write information about problem to file
  cout << "   ***Writing to file..." << endl;
  std::string filename_forward("Acoustic_HK_forward_LAR_N_YZ");
  Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename_forward);
  Write_Source(Proc_Idx, IW, IncWave::Pt, filename_forward);
  Write_Location(Proc_Idx, ScL, filename_forward);
  Write_Solution( Proc_Idx, Index, u, filename_forward);
  cout << "   ***Writing to file: done" << endl << endl;
  




 
  // *** Time-reversed field

  // Evaluate field on spherical manifold in far field

  std::vector<double> theta(size_t);
  std::vector<double> phi(size_p);
  std::vector<double> w_t(size_t);
  std::vector<double> w_p(size_p);
  getGaussLegendreQuad(size_t,theta,w_t);
  double val = 0.;
  for( int i_t = 0; i_t < size_t; i_t++ ){
    theta[i_t] = acos(theta[i_t]);
    //w_t[i_t] *= R*R;
  }
  //cout << "integration val : " << val << " : " << std::setprecision(17) << (val - (2.)/(81.)) << endl;
  
  for( int i_p = 0; i_p < size_p; i_p++ ){
    phi[i_p] = 2*PI*i_p/size_p;
    w_p[i_p] = 2*PI/size_p;
  }



  cout << " Preparing time-reversal migration..." << endl;
  std::vector<complex> src_rv(size_t*size_p);
  int k = 0;
  for( int i_t = 0; i_t < size_t; i_t++ ){
    cout << "i_t : " << i_t << endl;
    for( int i_p = 0; i_p < size_p; i_p++ ){
      Pvec p(R, theta[i_t], phi[i_p], Pvec::Spherical);

      // Evaluate expansions
      src_rv[k] = w_t[i_t]*w_p[i_p] * std::conj(Scatterer::Evaluate(p, ScL, u, Hankel,K_OUT) + IW->Evaluate(p));
      //src_rv[k] = w_t[i_t]*w_p[i_p] * std::conj( IW->Evaluate(p) );
      //cout << "src_rv : " <<  std::conj( IW->Evaluate(p)) << endl;

      k++;
    }
  }

  // Construct sources on spherical manifold in far field
  std::vector<IncWave*> IW_rv(size_t*size_p); 
  k = 0;
  std::vector<complex> expansion(Index.size());
  for( int i_t = 0; i_t < size_t; i_t++ ){
    for( int i_p = 0; i_p < size_p; i_p++ ){
      Pvec location(R,theta[i_t],phi[i_p],Pvec::Spherical);
      IW_rv[k] = new SphericalWave(K_OUT, src_rv[k], location);

      k++;
    }
  }
 
  

  std::vector< std::vector<complex> > RHS_rv(NScat, std::vector<complex>(Index.size()));
  std::vector< std::vector<complex> > u_rv(NScat, std::vector<complex>(Index.size()));
  cout << " Preparing time-reversal migration: done" << endl;
  
  
  
  cout << " Computing solution time-reversal problem..." << endl;
  // Apply H-K identity for Pressure waves
  res = 1e10; rel_res = 1e10;
  solver.RHSInit(IW_rv, ScL, RHS);  
  solver.Solve(PSC_env, ScL, RHS, u_rv, res, rel_res, cond);
  cout << "solved" << endl;

  // Write information about problem to file
  cout << "   ***Writing to file..." << endl;
  std::string filename_reverse("Acoustic_HK_reverse_LAR_N_YZ");
  Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename_reverse);
  Write_Source(Proc_Idx, IW_rv, IncWave::Pt, filename_reverse);
  Write_Location(Proc_Idx, ScL, filename_reverse);
  Write_Solution( Proc_Idx, Index, u_rv, filename_reverse);
  cout << "   ***Writing to file: done" << endl;


  // Clear memory
  //PSC_env.Destruct();



  // Post-processing (plot image and compute error)
  
  // TODO: Right now only x-y plane. Okay?

  cout << "   ***Creating images and computing error..." << endl;
  // Create a 2D grid onto which to compute field

  std::vector<double> X(2*M+1);
  for( int i = 0; i <= 2*M; i++ )
    X[i] = L * (i-M)/ ((double) M);

  // Compute image and store
  complex v_sc;
  complex v_rv;
  complex v_sc_in;
  complex v_rv_in;

  std::string str("Acoustic_Im_forward_Test_LAR_N_YZ.csv");
  char *Filename = (char*) str.c_str();
  ofstream im_f(Filename, ios::out);
 
  str = "Acoustic_Im_reverse_Test_LAR_N_YZ.csv";
  Filename = (char*) str.c_str();
  ofstream im_rv(Filename, ios::out);

  str = "grid.csv";
  Filename = (char*) str.c_str();
  ofstream grid(Filename, ios::out);



  // TODO: Transfer expansion from all transducers to the origin to obtain single expansion
  double err = 0;
  double norm = 0;
  for( int i = 0; i < (int) X.size(); i++ ){
    //cout << "row : " << i << endl;
    for( int j = 0; j < (int) X.size(); j++ ){
      //cout << "col : " << j << endl;
      Pvec p(X[i],X[j],0.,Pvec::Cartesian);

      bool too_close = false;
      Pvec q,s;
      for( int n = 0; n < NScat; n++){
	q = ScL[n].getLocation();
	s = p - q;
	if( s.r < 1.1*RADIUS ){
	  too_close = true;
	  break;
	}
      }

      // Record location of gridpoints
      grid << X[i] << "," << X[j] << "," << 0. << endl;
      
      // Evaluate scattered field at point

      // **** Forward field ****

      // Forward field P
      v_sc = 0;
      v_sc_in = 0.;
      if( !too_close ){
	v_sc_in = IW->Evaluate(p);
	v_sc = Scatterer::Evaluate(p, ScL, u, Hankel, K_OUT);
      }
      im_f << std::setprecision(15) << std::real(v_sc + v_sc_in) << "," << std::setprecision(15) << std::imag(v_sc + v_sc_in) << endl;
      //im_f << std::setprecision(15) << std::real( v_sc_in) << "," << std::setprecision(15) << std::imag( v_sc_in) << endl;

      // **** Time-reversed field P ****
      v_rv = 0;
      v_rv_in = 0.;
      if( !too_close ){
	for( int k = 0; k < (int) IW_rv.size(); k++){
	  complex v = IW_rv[k]->Evaluate(p);
	  v_rv_in += v;
	}

	v_rv = Scatterer::Evaluate(p, ScL, u_rv, Hankel, K_OUT);
      }
      //im_rv << std::setprecision(15) << std::real(v_rv_in) << "," << std::setprecision(15) << std::imag(v_rv_in) << endl;
      im_rv << std::setprecision(15) << std::real(v_rv + v_rv_in) << "," << std::setprecision(15) << std::imag(v_rv + v_rv_in) << endl;
      

      
    }
  }

  // Close files
  im_f.close();
  im_rv.close();
  grid.close();
  cout << "   ***Creating images and computing error: done" << endl;  




  





  return 0;
}






#endif
