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
//#include "./FMMPS/TransferUtils.h"

#include "./Aux/Eigen/Dense"


int main(int argc,char **args)
{

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
  std::string filename("../AcousticOutputFiles/FFH/HomogenizationCheck_O_");
  std::string response_filename("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_O_");

  std::string filename_FL("../AcousticOutputFiles/FFH/HomogenizationCheck_FL_O_");
  std::string response_filename_FL("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_FL_O_");

  std::string filename_single("../AcousticOutputFiles/FFH/HomogenizationCheck_single_O_");
  std::string response_filename_single("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_single_O_");



  // *** Parameters *** 
  double S = 1.; // Omega (cube) half-sidelength (Choosen such that Vol(Omega) = 1)
  double t = 1./3.;
  double delta = 0.75;    // For Eq.(1.8) 
  double R = 0.5;       // Cluster radius
  int K = 4;              // Number of satellites
  int J = 25;
  int Z = 30;            // Number of independent samples (for averaging)
  

  // *** Source parameters ***
  std::string src_type("PlaneWave");
  Pvec src_wv(0., 0., std::abs(K_OUT), Pvec::Cartesian);   // Plane wave, eave vector


  // *** Sampling points ***
  int nlevels = 3;
  int L = 15;             // Max degree expansion for translation at origin
  Indexing Index_FF(L);
  int N_S_MAX = 30;
  int N_S = 0;
  std::vector<Pvec> samples;

  for( int i = 0; i < N_S_MAX; i++ )
    samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, 0, Pvec::Spherical) );
  
  //for( int i = 0; i < N_S_MAX; i++ )
  // samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, PI/2., Pvec::Spherical) );
  
  N_S = samples.size();
 



  double a = 1e-2;
  for( int j = 0; j < J; j++){

    std::stringstream j_idx;
    j_idx << j; 
    
    a /= 1.5;
    double r = delta*a/2. / (K+1.);
    double V = 12.*a;
    int M = std::ceil( 1./V );          // Number of sub-domains Omega_m
    cout << "M : " << M << endl;
    //double s = (4*PI*pow(R, 3.)/3. / M) / 2;  // Omega_m (cube) half-sidelength (volume Omega_m = ALPHA*a)
    double s = pow(V, 1./3.)/2.;  // Omega_m (cube) half-sidelength (volume Omega_m = ALPHA*a)
    double d_min = std::pow(a, t);      // Minimum distance between 2 scatterers
    
    cout << "r :" << r << endl;
    cout << "s :" << s << endl;
    cout << "d_min :" << d_min << endl;
    
    
    complex C_in =  2*PI / pow(K_OUT*K_OUT - (K+1.)*4*PI*r/V, 1./2.) ;
    cout << "C_in : " << C_in << endl;
    
    
     // *** Initialization of source ***
      PlaneWave* IW = new PlaneWave(1., src_wv);
      
      
      
      for( int z = 0; z < Z; z++ ){

	std::stringstream z_idx;
	z_idx << z;
	
	// *** Construct scatterers ***
	cout << "***Building scatterers..." << endl;
	
	// Sets Omega_m (disjoint cubes with center gamma and sidelength c)
	srand48(time(0));
	int N = std::ceil( pow(M,1./3.)/2. );
	cout << "N : " << N << endl;
	M = (2*N)*(2*N)*(2*N);
	if( M > 5000 ) { nlevels = 4; }
	if( M > 20000 ) { nlevels = 5; }
	std::vector<Pvec> O_m(M);          // Centers of cubes Omega_m
	std::vector<Pvec> p_m(M);          // Location of "principal" scatterer in Omega_m
	std::vector< std::vector<Pvec> > s_m(M, std::vector<Pvec>() );          // Location of "scattelite" scatterers in Omega_m
	// *** Assumption: [K(z_m) + 1] = K(z_m) + 1 ***
	std::vector<int> K_m(M);      // Number of scatterers (K(m) + 1) in each sub-domains Omega_m
	std::vector<complex> C(M);
	int l = 0;
	std::vector<Scatterer> ScL;
	std::vector<Scatterer> ScL_FL;
	for( int i = -N; i < N; i++ ){
	  for( int j = -N; j < N; j++ ){
	    for( int k = -N; k < N; k++ ){
	      
	      K_m[l] = K;     // Number of satellites in Omega_m
	      
	      // TODO: Explain where does this 2i factor comes from. (Difference between green function 
	      // and Hankel function). minus sign comes from conjugate in formula in paper
	      C[l] = -CI/2.*C_OUT * (4*PI*r);  // Capacitance of scatterers (same for all) in Omega_m
	      
	      // TODO: Adjust so fast enough from side
	      O_m[l] = Pvec( i*(2.*s) + s, j*(2.*s) + s, k*(2.*s) + s, Pvec::Cartesian );
	      //std::vector<Pvec> v = RandCubic(s-r-d_min/2, s-r-d_min/2, s-r-d_min/2, r, d_min, O_m[l], K_m[l] + 1);
	      std::vector<Pvec> v(2);
	      v[0] = Pvec(O_m[l].x + d_min/2., O_m[l].y, O_m[l].z, Pvec::Cartesian);
	      v[1] = Pvec(O_m[l].x - d_min/2., O_m[l].y, O_m[l].z, Pvec::Cartesian);
	      
	      
	      //cout << i << " / " << N << endl;
	      // Full system (Spherical cluster)
	      if( O_m[l].r <= R ){
		p_m[l] =  v[0];
		Scatterer scatterer(r, K, K_OUT, RHO, p_m[l], 2);
		ScL.push_back(scatterer);
		
		// Foldy-Lax scatterers
		//Scatterer scatterer_FL(r, K, K_OUT, RHO, p_m[l], 2);
		//scatterer_FL.TM[0] = C[l] ;
		//for( int i = 1; i <= LMAX; i++ )
		//  scatterer_FL.TM[i] = 0.;
		//ScL_FL.push_back(scatterer_FL);
		
		
		s_m[l].resize(K_m[l]);
		for( int z = 0; z < K_m[l]; z++ ){ 
		  s_m[l][z] =  v[z+1];
		  Scatterer scatterer(r, K, K_OUT, RHO, s_m[l][z], 2);
		  ScL.push_back(scatterer);
		  
		  // Foldy-Lax scatterers
		  //Scatterer scatterer_FL(r, K, K_OUT, RHO, s_m[l][z], 2);
		  //scatterer_FL.TM[0] = C[l] ;
		  //for( int i = 1; i <= LMAX; i++ )
		  //  scatterer_FL.TM[i] = 0.;
		  //ScL_FL.push_back(scatterer_FL);
		}
	      }
	      
	      
	      l++;
	    }
	  }
	}
	
	
	
      
      cout << "***Building scatterers: done" << endl << endl;
      cout << "Total number of scatterers : " << ScL.size() << endl;
      
      
      
      
      
      
      // ---- Computations ----- //
      
      
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
      std::string temp_filename =  response_filename;
      temp_filename.append(j_idx.str());
      temp_filename.append("_");
      temp_filename.append(z_idx.str());
      temp_filename.append(".csv");
      char *Filename = (char*) temp_filename.c_str();
      ofstream Resfile(Filename, ios::out);
      for( int i = 0; i < (int) samples.size(); i++ )
	Resfile << samples[i].theta << ", " << samples[i].phi << ", " 
		<< std::real(response[i]) << ", " << std::imag(response[i]) << endl;
      Resfile.close();
      cout << " done" << endl;

    }


      
      // Destroy environment
      PSC_env.Destruct();
      
    
    
    /*

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
    
    */



    
    
    
    // --**-- Equivalent system (single sphere) --**-- 


      cout << endl <<  "   ***Single-sphere response...";
    // Single sphere's scattered expansion
    Pvec location(0.,0.,0., Pvec::Cartesian);
    Scatterer scatterer(R, 2.*PI/C_in, K_OUT, 1., location, 0);
    std::vector<complex> coeff((L+1)*(L+1));
    for( int i = 0; i < (L+1)*(L+1); i++ ){
      complex t_m = Scatterer::Mie(2*PI/C_in, K_OUT, R, 1, 0, Index_FF(i,0), Index_FF(i,1));
      coeff[i] = t_m * 4.*PI *  pow(CI,(double) Index_FF(i,0)) *
        std::conj( gsl_sf_harmonic(Index_FF(i,0), Index_FF(i,1),  src_wv.theta,  src_wv.phi) );
    }
    
    
    
    std::vector<complex> response_S(samples.size());
    for( int s = 0; s < (int) samples.size(); s++ ){
      for( int i = 0; i < ((L+1)*(L+1)); i++ ){
	int l = Index_FF(i,0);
	int m = Index_FF(i,1);
	
	response_S[s] += 1./(K_OUT*pow(CI, l+0.)) * coeff[i] * gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
      }
    }
    cout << "done" << endl;
    
    
    // Write response
    cout << "   ***Writing response...";
    std::string temp_filename_single =  response_filename_single;
    temp_filename_single.append(j_idx.str());
    temp_filename_single.append(".csv");
    char *Filename3 = (char*) temp_filename_single.c_str();
    ofstream Resfile_single(Filename3, ios::out);
    for( int i = 0; i < (int) samples.size(); i++ ){
      Resfile_single << samples[i].theta << ", " << samples[i].phi << ", " 
	      << std::real(response_S[i]) << ", " << std::imag(response_S[i]) << endl;
    }
    cout << " done" << endl;
    
  
    
    
    
    


  }
  
  return 0;
}





#endif
