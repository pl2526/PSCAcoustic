#ifndef MAIN_CPP
#define MAIN_CPP

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
#include "Imaging.h"
#include "Solver.h"
#include "Finalize.h"



int main(int argc,char **args)
{
  int Proc_Idx = 0;

  double begin = clock();
  for( int i = 0; i < 1.4e6; i++ ){
  }
cout << "Lowest Level Multipoles Done : " <<  (clock() - begin)*1./CLOCKS_PER_SEC << endl;

  // PSC environment for fast algorithm
  PSC_Env PSC_env(nLevels, EPS); 

  // Linear solver
  Solver solver; 

  // Global indexing
  Indexing Index(LMAX); // TODO: Should not have different instances of Index


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


  // *** Parameters ***

  int K = 50;  // Number of simulations to carry out

  // *** File names ***
  std::string filename("SuperRes_N=100000_A");
  std::string image_filename("SuperRes_N=100000_image_A");

  // *** Image parameters ***
  bool produce_image = false;
  int M = 70;     // Number of pixels: 2*M+1 from -M to M
  double L = 0.75;   // size of image [-L,L]x[-L,L]

  // *** Source parameters ***
  std::string src_type("PlaneWave");
  Pvec src_wv(0., 0., std::abs(K_OUT)*1., Pvec::Cartesian);   // Plane wave, eave vector

  // *** Scatterer cluster parameters *** 
  Pvec center(0.,0.,0.,Pvec::Cartesian);
  double d_min = 0.5*RADIUS;
  double R = 0.3;





  // ---- Computations -----

  double kappa = RHO / RHO_OUT;
  for( int k = 0; k < K; k++ ){
    Proc_Idx = k;

    std::vector< std::vector<complex > > RHS(NScat, std::vector<complex >(Index.size())); 
    std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
    std::vector< std::vector<complex> > u_in(NScat, std::vector<complex>(Index.size()));
    
    
    // *** Construct scatterers ***
    cout << "***Building scatterers..." << endl;
    std::vector<Scatterer> ScL;
    std::vector<Pvec> scat_loc = RandSpherical(R, RADIUS, d_min, center, NScat);
    assert(scat_loc.size() == NScat);
    for( int i = 0; i < NScat; i++ ){
      Scatterer scatterer(RADIUS, K, K_OUT, RHO, scat_loc[i]);
      ScL.push_back(scatterer);
    }
    cout << "***Building scatterers: done" << endl << endl;
    
    
    
    
    // *** Construct PSC environment ***
    cout << "***Constructing environment..." << endl;
    PSC_env.Construct(K_OUT, ScL);
    cout << "***Constructing environment: done" << endl << endl;
    
    
    
    // *** Initialization of source ***
    PlaneWave* IW = new PlaneWave(1., src_wv);
    
    
    
    // *** Initialization of right-hand side *** 
    cout << "***Initializing right-hand side..." << endl;
    solver.RHSInit(IW, ScL, RHS);  
    cout << "***Initializing right-hand side: done" << endl;
    
    
    
    // *** Solve linear system *** 
    cout << "   ***Solving linear system..." << endl;
    double res = 1e10, rel_res = 1e10, cond;
    
    solver.Solve(PSC_env, ScL, RHS, u, res, rel_res, cond);
    cout << "   ***Solving linear system: done" << endl;
    

    // *** Compute field inside scatterers *** 
    PSC_env.getTransfer()->execute(u, u_in);
    for( int n = 0; n < NScat; n++ ){
      std::vector<complex> u_src(Index.size());
      Pvec p = ScL[n].getLocation();
      IW->Translate(p, u_src);
      for( int i = 0; i < (int) Index.size(); i++ )
	u_in[n][i] += u_src[i];
    }

    for( int n = 0; n < NScat; n++ ){
      for( int i = 0; i < (int) Index.size(); i++ ){
	int l = Index(i,0);
	complex T_coeff = -K_OUT*RADIUS*RADIUS*CI * ( K_OUT*Amos_sf_hankel_l_prime(l, K_OUT*RADIUS) * Amos_sf_bessel_jl(l, K*RADIUS)
					 - K/kappa * Amos_sf_hankel_1(l, K_OUT*RADIUS) * Amos_sf_bessel_jl_prime(l, K*RADIUS) );
	T_coeff = pow(T_coeff, -1.);
	u_in[n][i] *= T_coeff;
      }
    }
    


    // *** Write information about problem to file ***
    cout << "   ***Writing to file..." << endl;
    Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename);
    Write_Source(Proc_Idx, IW, IncWave::Pt, filename);
    Write_Location(Proc_Idx, ScL, filename);
    Write_Solution( Proc_Idx, Index, u, filename);
    cout << "   ***Writing to file: done" << endl << endl;
    
    
    
    // *** Produce image ***
    std::vector< std::vector<complex> >* u_in_ptr; u_in_ptr = &u_in;
    Imaging(ScL, IW, image_filename, u, Proc_Idx, M, L, u_in_ptr);
    

    // Destroy environment
    PSC_env.Destruct();
  }




  return 0;
}





#endif
