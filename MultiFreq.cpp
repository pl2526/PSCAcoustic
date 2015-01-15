#ifndef MULTIFREQ_CPP
#define MULTIFREQ_CPP

/*
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved. 
 *
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

  // MPI Initialization
  int rank, size;
  MPI_Init(&argc, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &size ); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );


  // PSC environment for fast algorithm
  PSC_Env PSC_env(EPS); 

  // Linear solver
  Solver solver; 

  // Global indexing
  Indexing Index(LMAX); // TODO: Should not have different instances of Index



  // *** File names ***
  std::string filename("../AcousticOutputFiles/SR/SuperRes_B_PHI=0,05_N=2000_");
  std::string image_filename("../AcousticOutputFiles/SR/SuperRes_B_PHI=0,05_N=2000_image_");

  // *** Image parameters ***
  bool produce_image = true;
  int M = 65;     // Number of pixels: 2*M+1 from -M to M
  double L = 0.5;   // size of image [-L,L]x[-L,L]

  // *** Source parameters ***
  std::string src_type("PointSource");
  Pvec src_loc(0.05, 0.03, 0.,Pvec::Cartesian);   // Point source

  // *** Scatterer cluster parameters *** 
  Pvec center(0.,0.,0.,Pvec::Cartesian);
  double d_min = 0.002 * 2*PI/std::abs(K_OUT);
  double R = 0.3;
  double t = 0.2;
  double thick = 0.05;




  // *** Initialization *** //


  // Source signature
  double T = 80;                          // Imaging period
  double SIGMA = 1.;                     // Gaussian source
  double PHI = 1. / (2.*T);               // Sampling rate in frequency
  double CENTRAL = 1.;
  double a = 1;
  double S = 1./(a * 1./SIGMA);         // Sampling rate in time
  double F_MIN = -0.7;//-a * 1./SIGMA;                    // Min source frequency content
  F_MIN = std::max(PHI, std::floor( (CENTRAL + F_MIN)/PHI ) * PHI);    // (so it is a multiple of S);
  double F_MAX = 0.7;//a * 1./SIGMA;           // Max source frequency content
  F_MAX = std::ceil( (CENTRAL + F_MAX)/PHI ) * PHI;     // (so it is a multiple of S);

  
  int F_local_size = std::ceil( ( (double) F_MAX - (double) F_MIN ) / PHI / (double) size);
  cout << "local size : " << F_local_size << endl;
  int F_global_size = size*F_local_size;
  double* freq_global = new double[ F_global_size ]; 
  
  if( rank == 0 ){

  
    std::string extra_filename = filename + "extra.csv";
    ofstream Extrafile((char*) extra_filename.c_str(), ios::out);
    
    // Index
    Extrafile << " Period [-T;T] :," << T << endl;
    Extrafile << " Standard deviation (sigma) :," << SIGMA << endl;
    Extrafile << " Sampling in frequency :," << PHI << endl;
    Extrafile << " Sampling in time :," << S << endl;
    Extrafile << " Central frequency :," << CENTRAL << endl;
    Extrafile << " Minimum nonzero frequency :," << F_MIN << endl;
    Extrafile << " Maximum nonzero frequency :," << F_MAX << endl;
    Extrafile << " Number of pixels in each dimension :," << 2*M+1 << endl;
    
    Extrafile.close();

    for( int i = 0; i < F_global_size; i++ ){
      double freq = F_MIN + i*PHI;
      freq_global[i] = freq;
    } 

  // Display relevant information concerning the problem
    cout << "Number of processes : " << size << endl;
    cout << "Number of tasks : " << F_global_size << endl;
    cout << "F_MIN : " << F_MIN << endl;
    cout << "F_MAX : " << F_MAX << endl;
    cout << endl << endl;
    cout << "LMAX : " << LMAX << endl;
    cout << "RTOL : " << RTOL << endl;
    cout << "MAXITS : " << MAXITS << endl;
    cout << "OMEGA : " << OMEGA << endl;
    cout << "RHO_OUT : " << RHO_OUT << endl;
    cout << "RHO : " << RHO << endl;
    cout << "C_OUT : " << C_OUT << endl;
    cout << "C : " << C << endl;
    cout << "K_OUT : " << K_OUT << endl;
    cout << "K : " << K << endl;
    cout << "NScat : " << NScat << endl;
    cout << "RADIUS : " << RADIUS << endl;
    cout << "nLevels : " << nLevels << endl;
    cout << "EPS : " << EPS << endl;
    cout << endl << endl;

  }
    


  // ---- Parallelization ----
  // Distribute frequencies
  double* freq_local = new double[ F_local_size ];

  // Distribute problems for different frequencies among processes
  MPI_Scatter(&freq_global[0], F_local_size, MPI_DOUBLE, &freq_local[0] , F_local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for( int i = 0; i < F_local_size; i++ )
    cout << size << " : " << rank << " : " << freq_local[i] << endl;
  MPI_Barrier( MPI_COMM_WORLD );

  

 // ---- Computations -----
  double res = 1e10, rel_res = 1e10, cond;
  int Niter;
  for( int n = 0; n <  F_local_size; n++ ){
    int proc_idx = std::floor(freq_local[n] / PHI + 0.5);

    complex k_in = 2.*PI* freq_local[n] / C;
    complex k_out = 2.*PI* freq_local[n] / C_OUT;

    cout << rank << ":" << freq_local[n] << " : " << k_in << " : " << k_out << endl;


    std::vector< std::vector<complex > > RHS(NScat, std::vector<complex >(Index.size())); 
    std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
    std::vector< std::vector<complex> > u_in(NScat, std::vector<complex>(Index.size()));
    
    // !!!!! NEED TO MAKE SURE RANDOM SEED IS ALL SAME
    // *** Construct scatterers ***
    cout << "***Building scatterers..." << endl;
    std::vector<Scatterer> ScL;
    std::vector<Pvec> scat_loc = ErgodicCavity(R, t, thick, center, RADIUS, d_min, NScat, src_loc, false);
    //std::vector<Pvec> scat_loc = TruncatedSphere(R, t, RADIUS, d_min, center, NScat, src_loc, false);
    assert(scat_loc.size() == NScat);
    for( int i = 0; i < NScat; i++ ){
      Scatterer scatterer(RADIUS, k_in, k_out, RHO, scat_loc[i]);
      ScL.push_back(scatterer);
    }
    cout << "***Building scatterers: done" << endl << endl;

    // *** Initialization of source ***
    cout << "***Constructing source..." << endl;
    SphericalWave* IW = new SphericalWave(k_out, 1., src_loc);
    cout << "***Constructing source: done" << endl << endl;
    
    
      
    
    // *** Construct PSC environment ***
    cout << "***Constructing environment..." << endl;
    PSC_env.Construct(nLevels, RADIUS, k_in, k_out, ScL);
    cout << "***Constructing environment: done" << endl << endl;
    
    
    
    
    
    // *** Initialization of right-hand side *** 
    cout << "***Initializing right-hand side..." << endl;
    solver.RHSInit(IW, ScL, RHS);  
    cout << "***Initializing right-hand side: done" << endl;
    
    
    
    // *** Solve linear system *** 
    cout << "   ***Solving linear system..." << endl;
    solver.Solve(PSC_env, ScL, RHS, u, res, rel_res, cond, Niter, proc_idx);
    cout << "   ***Solving linear system: done" << endl;
    

    // *** Compute expansion inside scatterers *** 
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
	u_in[n][i] /= Scatterer::Mie_in(k_in, k_out, RADIUS, RHO, RHO_OUT, l);
      }
   }
    


    // *** Write information about problem to file ***
    cout << "   ***Writing to file..." << endl;
    Write_Info(proc_idx, src_type, res, rel_res, cond, Niter, filename, k_out, k_in);
    Write_Source(proc_idx, IW, IncWave::Pt, filename);
    Write_Location(proc_idx, ScL, filename);
    Write_Solution( proc_idx, Index, u, filename);
    cout << "   ***Writing to file: done" << endl << endl;
    
    // *** Produce image ***
    if( produce_image ){
      std::vector< std::vector<complex> >* u_in_ptr; u_in_ptr = &u_in;
      Imaging(ScL, IW, image_filename, u, proc_idx, M, L, u_in_ptr, k_out, k_in);
    }
    

    // Destroy environment
    delete IW;
    PSC_env.Destruct();
  }

  // Wait before finalizing
  MPI_Barrier( MPI_COMM_WORLD );
  
  // Clear memory
  MPI_Finalize();

  return 0;
}





#endif
