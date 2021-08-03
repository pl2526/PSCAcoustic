#ifndef MAIN_CPP
#define MAIN_CPP

/*
 *  MAIN.h
 *  PSC
 *
 *
 *  Overstructure.
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


  // !!!!! Make sure Mie coefficients are those of sound-soft sphere


  // *** Parameters ***



  // *** File names ***
  std::string filename("../AcousticOutputFiles/HomogenizationCheck_F_PHI=0,12_N=5000_");
  std::string image_filename("../AcousticOutputFiles/HomogenizationCheck_F_PHI=0,12_N=5000_image_");


  // *** Source parameters ***
  std::string src_type("PlaneWave");
  Pvec src_wv(0., std::abs(K_OUT)*1., 0., Pvec::Cartesian);   // Plane wave, eave vector

  // *** Scatterer cluster parameters *** 
  Pvec center(0.,0.,0.,Pvec::Cartesian);
  double R = 0.4;

  // *** Sampling points ***
  int N_S = 1;
  std::vector<Pvec> samples(N_S);
  samples(1) = Pvec(1., PI/2., PI/2., Pvec::Spherical);
  



  // ---- Computations -----
  
  Proc_Idx = 0;
  
  std::vector< std::vector<complex > > RHS(NScat, std::vector<complex >(Index.size())); 
  std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
  
  // *** Construct scatterers ***
  cout << "***Building scatterers..." << endl;

  double S = 1./2.; // Omega (cube) half-sidelength (Choosen such that Vol(Omega) = 1)
  double t = 1./3.;
  doube a = 0.1;
  int M = std::ceil( 1./a );          // Number of sub-domains Omega_m
  double s = pow(a, 1./3.)/2.;  // Omega_m (cube) half-sidelength (volume Omega_m = a)
  double d_min = std::pow(a, t);      // Minimum distance between 2 scatterers

  
  // *** Assumption: [K(z_m) + 1] = K(z_m) + 1 ***
  std::vector<int> K(M);      // Number of scatterers (K(m) + 1) in each sub-domains Omega_m
  std::vector<double> sigma(M);
  for( int m = 0; m < M; m++ ){
    K[m] = 0;
    sigma[m] = 1.;
  }

  
  // Sets Omega_m (disjoint cubes with center gamma and sidelength c)
  int N = std::ceil(S / (2*s) );
  std::vector<Pvec> O_m(M);
  int l = 0;
  for( int i = -N; i < N; i++ ){
    for( int j = -N; j < N; j++ ){
      for( int k = -N; k < N; k++ ){
	O_m[l] = Pvec( i*(2.*s) + s, j*(2.*s) + s, k*(2.*s) + s, Pvec::Cartesian );
	cout << O_m[l].x << " : " <<  O_m[l].y << " : " <<  O_m[l].z << endl;
	l++;
      }
    }
  }
  cout << "l : " << M << endl
 
  




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
  
  
  
  
  // *** Write information about problem to file ***
  cout << "   ***Writing to file..." << endl;
  Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename);
  Write_Source(Proc_Idx, IW, IncWave::Pt, filename);
  Write_Location(Proc_Idx, ScL, filename);
  Write_Solution( Proc_Idx, Index, u, filename);
  cout << "   ***Writing to file: done" << endl << endl;
  
  
  

  // *** Compute far-field signature ***

  // Translate to origin
  int L = 10;
  std::vector< std::vector<complex> > response(N_S, std::vector<complex>((L+1)*(L+1)) );
  std::vector<complex> FF((L+1)*(L+1));
  Indexing Index_FF(L);
  for( int n = 0; n < NScat; n++ ){
    Pvec loc = -ScL[n].getLocation();

    FMM::PointAndShoot T(loc.r, K_OUT, loc.theta, loc.phi, Index_FF, Index, true);
    std::vector<complex> vec = T.Apply(u[n]);

    for( int i = 0; i < (L+1)*(L+1); i++)
      FF[i] += vec[i];
  }

  // TODO: Initialize reponse?
  for( int s = 0; samples.size(); s++ ){
    for( int i = 0; i < (L+1)*(L+1); i++ ){
      int l = Index_FF(i,0);
      int m = Index_FF(i,1);
      response[s][i] += 1./(K_OUT*pow(CI, l+0.)) * FF[i] * gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
    }
  }
  
  // Destroy environment
  PSC_env.Destruct();





  return 0;
}





#endif
