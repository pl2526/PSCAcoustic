#ifndef MAIN_CPP
#define MAIN_CPP

/*
 *  MAIN.h
 *  PSC
 *
 *  Overstructure.
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
//#include "./FMMPS/TransferUtils.h"

#include "./Aux/Eigen/Dense"


int main(int argc,char **args)
{
  int Proc_Idx = 0;

  // PSC environment for fast algorithm
  PSC_Env PSC_env(EPS); 

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



  // *** File names ***
  std::string filename("../AcousticOutputFiles/FFH/HomogenizationCheck_B_");
  std::string response_filename("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_B_");

  std::string filename_FL("../AcousticOutputFiles/FFH/HomogenizationCheck_FL_B_");
  std::string response_filename_FL("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_FL_B_");



  // *** Parameters *** 
  double s = 1.; // Omega (cube) half-sidelength (Choosen such that Vol(Omega) = 1)
  double t = 1./3.;
  double R = 0.5;       // Cluster radius

  

  // *** Source parameters ***
  std::string src_type("PlaneWave");
  Pvec src_wv(0., 0., std::abs(K_OUT), Pvec::Cartesian);   // Plane wave, eave vector


  // *** Sampling points ***
  int nlevels = 3;
  int L = 10;             // Max degree expansion for translation at origin
  Indexing Index_FF(L);
  int N_S_MAX = 30;
  int N_S = 0;
  std::vector<Pvec> samples;

  for( int i = 0; i < N_S_MAX; i++ )
    samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, 0, Pvec::Spherical) );
  
  //for( int i = 0; i < N_S_MAX; i++ )
  // samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, PI/2., Pvec::Spherical) );
  
  N_S = samples.size();
 

  int J = 40;
  double a = 1e-2;
  for( int j = 0; j < J; j++){

    a /= 1.5;
    double d_min = std::pow(a, t);      // Minimum distance between 2 scatterers
    int M = std::ceil( 1./pow(a, s) );          // Number of sub-domains Omega_m
    double V = 2.85 * M * 4*PI/3*pow(a+d_min/2,3.);
    R = pow(V,1./3.)/2.;
    cout << "Number of scatterers : " << M << endl;

    complex C_in =  2*PI / pow(K_OUT*K_OUT - 4*PI, 1./2.) ;
    cout << "C_in : " << C_in << endl;

    
    // *** Construct scatterers ***
    cout << "***Building scatterers..." << endl;
    
    // Sets Omega_m (disjoint cubes with center gamma and sidelength c)
    srand48(time(0));
    if( M > 5000 ) { nlevels = 4; }
    if( M > 50000 ) { nlevels = 5; }
    std::vector<Scatterer> ScL;
    std::vector<Scatterer> ScL_FL;
    Pvec loc( 0.,0.,0., Pvec::Cartesian );
    std::vector<Pvec> v = RandCubic(R, R, R, a, d_min, loc, M);
    for( int i = 0; i < (int) v.size(); i++ ){

      Scatterer scatterer(a, K, K_OUT, RHO, v[i], 2);
      ScL.push_back(scatterer);
      
      // Foldy-Lax scatterers
      // TODO: Explain where does this 2i factor comes from. (Difference between green function 
      // and Hankel function). minus sign comes from conjugate in formula in paper
      Scatterer scatterer_FL(a, K, K_OUT, RHO, v[i], 2);
      scatterer_FL.TM[0] = -2.*PI* a * CI;
      for( int i = 1; i <= LMAX; i++ )
	scatterer_FL.TM[i] = 0.;
      ScL_FL.push_back(scatterer_FL);
      
    }
	    


    
    cout << "***Building scatterers: done" << endl << endl;
    cout << "Total number of scatterers : " << ScL.size() << endl;
    
    



    
    // ---- Computations ----- //

    // *** Initialization of source ***
    PlaneWave* IW = new PlaneWave(1., src_wv);
    
    
    // --**-- Full system --**--

    std::vector< std::vector<complex > > RHS(ScL.size(), std::vector<complex >(Index.size())); 
    std::vector< std::vector<complex> > u(ScL.size(), std::vector<complex>(Index.size()));

    // *** Construct PSC environment ***
    cout << "***Constructing environment...";
    PSC_env.Construct(nlevels, a, K, K_OUT, ScL);
    cout << " done" << endl << endl;
    
    
    // *** Initialization of right-hand side *** 
    cout << "***Initializing right-hand side..." << endl;
    solver.RHSInit(IW, ScL, RHS);  
    cout << "***Initializing right-hand side: done" << endl;
    
    
    // *** Solve linear system *** 
    cout << "   ***Solving linear system..." << endl;
    double res = 1e10, rel_res = 1e10, cond;
    int Niter;
    
    solver.Solve(PSC_env, ScL, RHS, u, res, rel_res, cond, Niter);
    cout << "   ***Solving linear system: done" << endl;
    
    
    
    // TODO: Write in different format?
    // *** Write information about problem to file ***
    cout << "   ***Writing to file...";
    Write_Info(j, src_type, res, rel_res, cond, Niter, filename, K_OUT, K, ScL.size());
    Write_Source(j, IW, IncWave::Pt, filename);
    Write_Location(j, ScL, filename);
    Write_Solution(j, Index, u, filename);
    cout << " done" << endl << endl;
    
    
    
    
    // *** Compute far-field signature ***
    cout << "   ***Computing far-field..." << endl;
    // Translate to origin
    std::vector<complex> response(N_S);

    // TODO: Initialize reponse?
    double dot;
    for( int s = 0; s < (int) samples.size(); s++ ){
	for( int i = 0; i < ((LMAX+1)*(LMAX+1)); i++ ){
	  int l = Index(i,0);
	  int m = Index(i,1);
	  
	  
	  //response[s] += 1./(K_OUT*pow(CI, l+0.)) * FF[i] * gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
	  for( int n = 0; n < (int) ScL.size(); n++ ){
	    Pvec vec = -ScL[n].getLocation();
	    dot = vec.x*samples[s].x + vec.y*samples[s].y + vec.z*samples[s].z;
	    response[s] += 1./(K_OUT*pow(CI, l+0.)) * u[n][i] *  exp(CI*K_OUT*dot) *
	      gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
	  }

	}
    }
    cout << " done" << endl;
    
    
    cout << "   ***Writing response...";
    std::stringstream nProc_str;
    nProc_str << j;
    std::string temp_filename =  response_filename;
    temp_filename.append(nProc_str.str());
    temp_filename.append(".csv");
    char *Filename = (char*) temp_filename.c_str();
    ofstream Resfile(Filename, ios::out);
    for( int i = 0; i < (int) samples.size(); i++ ){
      Resfile << samples[i].theta << ", " << samples[i].phi << ", " 
	      << std::real(response[i]) << ", " << std::imag(response[i]) << endl;
      
      //cout << response[i] << endl;
    }
    Resfile.close();
    cout << " done" << endl;

    // Destroy environment
    PSC_env.Destruct();
    
    
    
    

    // --**-- Foldy-Lax system --**--
    
    cout << "***Constructing environment...";
    PSC_env.Construct(nlevels, a, K, K_OUT, ScL_FL);
    cout << " done" << endl << endl;


    // *** Initialization of right-hand side *** 
    cout << "***Initializing right-hand side...";
    solver.RHSInit(IW, ScL_FL, RHS);  
    cout << " done" << endl;
    
    
    // *** Solve linear system *** 
    cout << "   ***Solving linear system...";
    res = 1e10; rel_res = 1e10;
    solver.Solve(PSC_env, ScL_FL, RHS, u, res, rel_res, cond, Niter);
    cout << " done" << endl;
    
    
    
    // TODO: Write in different format?
    // *** Write information about problem to file ***
    cout << "   ***Writing to file...";
    Write_Info(j, src_type, res, rel_res, cond, Niter, filename_FL, K_OUT, K, ScL.size());
    Write_Source(j, IW, IncWave::Pt, filename_FL);
    Write_Location(j, ScL, filename_FL);
    Write_Solution( j, Index, u, filename_FL);
    cout << " done" << endl << endl;
    
    
        // *** Compute far-field signature ***
    cout << "   ***Computing far-field..." << endl;
    // Translate to origin
    std::vector<complex> response_FL(N_S);

    // TODO: Initialize reponse?
    for( int s = 0; s < (int) samples.size(); s++ ){
	for( int i = 0; i < ((LMAX+1)*(LMAX+1)); i++ ){
	  int l = Index(i,0);
	  int m = Index(i,1);

	  for( int n = 0; n < (int) ScL_FL.size(); n++ ){
	    Pvec vec = -ScL_FL[n].getLocation();
	    double dot = vec.x*samples[s].x + vec.y*samples[s].y + vec.z*samples[s].z;
	    response_FL[s] += 1./(K_OUT*pow(CI, l+0.)) * u[n][i] *  exp(CI*K_OUT*dot) *
	      gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
	  }

	}
    }
    cout << " done" << endl;

    
    
    cout << "   ***Writing response...";
    std::string temp_filename_FL =  response_filename_FL;
    temp_filename_FL.append(nProc_str.str());
    temp_filename_FL.append(".csv");
    char *Filename2 = (char*) temp_filename_FL.c_str();
    ofstream Resfile_FL(Filename2, ios::out);
    for( int i = 0; i < samples.size(); i++ ){
      Resfile_FL << samples[i].theta << ", " << samples[i].phi << ", " 
		    << std::real(response_FL[i]) << ", " << std::imag(response_FL[i]) << endl;
    }
    cout << " done" << endl;
    
        
    // Destroy environment
    PSC_env.Destruct();
    
    


  }
  
  return 0;
}





#endif
