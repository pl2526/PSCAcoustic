#ifndef DI_CPP
#define DI_CPP

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

#include "General.h"
#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "Imaging.h"
#include "PSCEnv.h"
#include "Solver.h"
#include "Finalize.h"

#include "./FMMPS/TransferUtils.h"

Indexing temp_index(LMAX);


// **** For single-scattering, don't forget to adjust the Solver.cpp
// **** Dont't forget to set L = 0

// Scattering coefficients of point-like satterers
double ScatCoeff(double mu, double a){
  double val;
  double x = drand48();
  if( x < 0.5 )
    val = mu + a;
  else
    val = mu - a;

  return val;
}

// Displacement profile
// *** Gradient of the displacement must be of the order of eps^(1+2*delta)
Pvec Displacement( Pvec p, double amp, double eps, double delta){
  double x_disp = 5e-4;// * cos(2*PI*p.x);//8e-3 * sin(PI*(p.x-2));//2e-3;//exp(-pow(eps, delta)/2. * p.x*p.x);
  double y_disp = 0;//2e-3 * sin(2*PI*p.y);//2e-3;//amp * 1.;
  double z_disp = 0.;
  Pvec disp(x_disp, y_disp, z_disp, Pvec::Cartesian);
  return disp;
}

// TODO: WHAT IS THIS LOG(1e-5)?
// Windows for curvelets
double W_r(double scale, double r){
    double center_r = 5./(4.*scale);
    double sigma_r = -pow(3./(4.*scale),2.) / log(1e-5);
    return std::exp(-pow(r - center_r,2.)/sigma_r);
}
   
double W_p(double scale, double center_p, double p){
  double sigma_p = -scale/log(1e-5);
  return std::exp(-pow(p - center_p,2.)/sigma_p);
}    


int main(int argc,char **args)
{
  int Proc_Idx = 0;

  //PSC_Env PSC_env(nLevels, EPS);  // PSC environment for fast algorithm
  Solver solver;   // Linear solver
  Indexing Index(LMAX); // TODO: Should not have different instances of Index
  Indexing local_Index(5);
  FMM::S_Rotation::S_Init(12);
  FMM::Z_transfer::Z_Init(12);

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


  // *** File names ***
  std::string curvelet_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_curvelet_H");

  std::string filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_H");
  std::string freq_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_freq_H");
  std::string extra_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_extra_H");
  std::string disp_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_H_disp");
  std::string src_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_src_H");
  std::string rcv_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_rcv_H");
  std::string sgn_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_sgn_H");
  std::string image_filename("../AcousticOutputFiles/DI/DI_SingleScattering_mu=0_image_H");

  std::string filename_d(filename); filename_d.append("_disp");
  std::string freq_filename_d(freq_filename); freq_filename_d.append("_disp.csv");
  std::string extra_filename_d(extra_filename); extra_filename_d.append("_disp.csv");
  std::string disp_filename_d(disp_filename); disp_filename_d.append("_disp.csv");
  std::string sgn_filename_d(sgn_filename); sgn_filename_d.append("_disp");
  std::string image_filename_d(image_filename); image_filename_d.append("_disp.csv");

  curvelet_filename.append(".csv");
  freq_filename.append(".csv");
  extra_filename.append(".csv");
  disp_filename.append(".csv");
  src_filename.append(".csv");
  rcv_filename.append(".csv");
  //sgn_filename.append(".csv");
  image_filename.append(".csv");


  // *** General parameters ***
  std::string src_type("PW");
  double mu = 0.;               // Mean of scattering coefficients
  double sigma = 1.;          // Variance of scattering coefficients
  double a = 1.;                // amplitude of scattering coefficients
  double src_scale = 5.e-3/(2*PI);      // Scale of the curvelet
  double delta = pow(src_scale, 1.5);
  double disp_amp_scale =  1e-4;//pow(src_scale, 1.+delta);           // Scale of the amplitude of the displacement
  int N_t = 200; // Number of samples in time



  // *** Sources' parameters ***
  int N_src = 1;
  std::vector<Pvec> src_dir(N_src);
  std::vector<Pvec> src_loc(N_src);
  src_dir[0] = Pvec(1.,PI/2.,0.,Pvec::Spherical);
  src_loc[0] = Pvec(0.,0.,0.,Pvec::Cartesian);

  // Save sources information
  ofstream src_file(src_filename.c_str(), ios::out);
  for( int q = 0; q < N_src; q++ ){
    src_file << src_dir[q].x << "," <<  src_dir[q].y << "," << src_dir[q].z << "," << src_loc[q].x << "," <<  src_loc[q].y << "," << src_loc[q].z << endl;
  }
  src_file.close();




// *** Receivers' locations ***
  int N_rcv = 100;
  double R = 1e6;
  std::vector<Pvec> rcv_loc(N_rcv);
  for( int i = 0; i < N_rcv; i++){
    double phi_rec = 2*PI * (double) i / (N_rcv + 1.);
    rcv_loc[i] = Pvec(5e3, PI/2., phi_rec, Pvec::Spherical);
  }
  /*std::vector<Pvec> rcv_loc(5);
  rcv_loc[0] =  Pvec(-5e3,0.,0.,Pvec::Cartesian);
  rcv_loc[1] =  Pvec(0.,-5e3,0.,Pvec::Cartesian);
  rcv_loc[2] =  Pvec(5e3/sqrt(2),5e3/sqrt(2), 0.,Pvec::Cartesian);
  rcv_loc[3] =  Pvec(5e3, -3e3,0.,Pvec::Cartesian);
  rcv_loc[4] =  Pvec(5e3, 3e3,0.,Pvec::Cartesian);*/

  // TODO: Needed?
  double max_rcv_loc_r = 0;
  for( int i = 0; i < rcv_loc.size(); i++ )
    max_rcv_loc_r = std::max(max_rcv_loc_r, rcv_loc[i].r);
  
  // Save receivers information
  ofstream rcv_file(rcv_filename.c_str(), ios::out);
  for( int q = 0; q < N_rcv; q++ ){
    rcv_file << rcv_loc[q].x << "," <<  rcv_loc[q].y << "," << rcv_loc[q].z << endl;
  }
  rcv_file.close();
 



  // *** Cluster properties ***
  double X_DIM = 1.1*sqrt(src_scale);
  double Y_DIM = 0.1;//0.7;
  Pvec Sc_center(0.2, 0., 0.,Pvec::Cartesian);
  double Sc_max = Sc_center.x + std::max(X_DIM, Y_DIM);
  

  

  // *** Time samples, curvelet location and displacement ***
  ofstream crv_file((char*) curvelet_filename.c_str(), ios::out);
  ofstream disp_file((char*) disp_filename.c_str(), ios::out);
  std::vector< std::vector<double> > time(N_rcv, std::vector<double>(N_t) );
  for( int i = 0; i < N_t; i++ ){

    // Location of a curvelet leaving from the origin in the x-direction
    Pvec loc( (double) i / (double) N_t * 1.2*Sc_max, 0., 0., Pvec::Cartesian);
    crv_file << loc.x << "," << loc.y << "," << loc.z << endl;

    // Displacement at location of curvelet
    Pvec disp = Displacement(loc, disp_amp_scale, src_scale, delta);
    disp_file << disp.x << "," << disp.y << "," << disp.z << endl;

    // Time at which a curvelet located at "loc" will arrive at receiver
    for(int q = 0; q < N_rcv; q++){
      Pvec x_hat = rcv_loc[q] - loc;
      time[q][i] = std::sqrt( (rcv_loc[q].x - loc.x)*(rcv_loc[q].x - loc.x) + (rcv_loc[q].y - loc.y)*(rcv_loc[q].y - loc.y) +
			      (rcv_loc[q].z - loc.z)*(rcv_loc[q].z - loc.z) ) + loc.r;
    }
  }
  crv_file.close();
  disp_file.close();


  



  // *** Frequency content of curvelet ***

   // xi_r corresponds to radius; xi_p corresponds to phi
  double src_rate = 1./(1. * Sc_max);  
  //double src_rate = 1./(6.*Y_DIM);  
    double xi_r_min = 1./(2.*src_scale);
    double xi_r_max = 2./src_scale;
    //double xi_p_min = - 1.1*sqrt(src_scale);
    //double xi_p_max =  1.1*sqrt(src_scale);
    //double xi_p_min = src_dir[s].phi - 1.1*sqrt(src_scale);
    //double xi_p_max = src_dir[s].phi + 1.1*sqrt(src_scale);
    
    // Sampling
    double x_min = -xi_r_max;
    double x_max = xi_r_max;
    int N_crv = 2*std::ceil( std::abs(x_max) / src_rate);
    cout << "Number of samples : " << N_crv << ":" << N_crv * N_crv << endl;
    
    double src_xi_x, src_xi_y;
    complex c;
    std::vector< std::vector<Pvec> > src_xi(N_src, std::vector<Pvec>() );
    std::vector<int> N_freq(N_src);
    std::vector<std::vector<std::vector<complex> > > C(N_src, std::vector<std::vector<complex> >(N_rcv, std::vector<complex>()));
    
    
    for( int s = 0; s < N_src; s++ ){
      int i = 0;
      int j = 0;
      int ctr = 0;
      
      while( i <= N_crv ){
	cout << i << endl;
	src_xi_x = x_min + ((double) i)/N_crv * (x_max - x_min);
	while( j <= N_crv ){
	  src_xi_y = x_min + ((double) j)/N_crv * (x_max - x_min);
	  
	  Pvec src_wv(src_xi_x, src_xi_y, 0., Pvec::Cartesian);	
	  c =  W_r(src_scale, src_wv.r) * W_p(src_scale, src_dir[s].phi, src_wv.phi);
	  
	  if( std::abs(c) > 1e-9 ){
	    Pvec vec(src_xi_x, src_xi_y, 0., Pvec::Cartesian);
	    src_xi[s].push_back(vec);
	    ctr++;
	    
	    for( int q = 0; q < N_rcv; q++){
	      C[s][q].push_back( exp(CI*src_wv.r*rcv_loc[q].r) * c );
	    }
	  }
	  
	  j++;
	}
	j = 0;
	i++;
      }
      
      N_freq[s] = ctr;
    }

    for( int s = 0; s < N_src; s++ )
      src_xi[s].resize(N_freq[s]);








  // *** Save extra information to file ***
    // TODO: FIX ME
  // Save extra info to file (info that's not used in the standard cases)
  char *Filename = (char*) extra_filename.c_str();
  ofstream extra_file(Filename, ios::out);
  extra_file << "Number of freq samples in x: ," << 0. << endl;
  extra_file << "Number of freq samples in y: ," << 0. << endl;
  extra_file << "Freq sampling rate (radius): ," << src_rate << endl;
  extra_file << "Freq sampling rate (phi): ," << src_rate << endl;
  extra_file << "Curvelet scale: ," << src_scale << endl;
  extra_file << "Curvelet direction: ," << src_dir[0].x << " , "<<  src_dir[0].y << " , " << src_dir[0].z << " , " << endl;
  extra_file << "Cluster center: ," << Sc_center.x << "," << Sc_center.y << "," << Sc_center.z  << endl;
  extra_file << "Cluster X dim: ," << X_DIM << endl;
  extra_file << "Cluster Y dim: ," << Y_DIM << endl;
  //extra_file << "Time interval: ," << T_1 << "," << T_2 << endl;
  extra_file << "Time interval samples: ," << N_t << endl;
  //extra_file << "Image window center: ," << N_d_x << endl;
  //extra_file << "Image window N_x pixels: ," << N_d_x << endl;
  //extra_file << "Image window N_y pixels: ," << N_d_y << endl;







  // ---- Computations -----
  // *** Compute solutions for each source and receiver locations ***

  std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
  double res, rel_res, cond;
  Proc_Idx = 0;
  std::vector< std::vector< std::vector<complex> > > trace(N_src,std::vector<std::vector<complex> >(N_rcv,std::vector <complex>(N_t,0.)));
 std::vector< std::vector< std::vector<complex> > > trace_disp(N_src,std::vector<std::vector<complex> >(N_rcv,std::vector <complex>(N_t,0.)));


  // Initialization
 // TODO: Needed?
 for( int s = 0; s < N_src; s++){
   for( int q = 0; q < N_rcv; q++ ){
     for( int j = 0; j < N_t; j++ ){
       trace[s][q][j] = 0.;
       trace_disp[s][q][j] = 0.;
     }
   }
 }


    
    
    // *** Construct scatterers ***
    cout << "***Building scatterers..." << endl;
    std::vector<Scatterer> ScL;
    std::vector<Scatterer> ScL_disp;
    double d_min = 0.;
    //std::vector<Pvec> scat_loc(1); scat_loc[0] = Pvec(0.,1.1,0.,Pvec::Cartesian);
    //std::vector<Pvec> scat_loc = RandSphericalXY(R, RADIUS, d_min, center, NScat);
    std::vector<Pvec> scat_loc = RandRectangularXY(X_DIM, Y_DIM, RADIUS, d_min, Sc_center, NScat);
    assert(scat_loc.size() == NScat);
    std::vector<complex> t_matrix_coeff(NScat); // Nonscales  
    for( int n = 0; n < NScat; n++ ){
      
      Pvec loc(scat_loc[n].x, scat_loc[n].y, 0., Pvec::Cartesian);
      Pvec loc_disp = loc + Displacement(loc, disp_amp_scale, src_scale, delta);
      cout << "loc : " <<  loc.x << " : " << loc.y << " : " << loc.z << " : " << endl;
      double s = 0;
      ( drand48() < 0.5 ) ? s = -sigma : s = sigma;
      t_matrix_coeff[n] = mu + s; 
      
      Scatterer scatterer( RADIUS, K, K_OUT, RHO, loc);
      ScL.push_back(scatterer);
      
      Scatterer scatterer_disp( RADIUS, K, K_OUT, RHO, loc_disp);
      ScL_disp.push_back(scatterer_disp);
      
    }
   
    for( int n = 0; n < NScat; n++ ){
      ScL[n].TM[0] = t_matrix_coeff[n];
      ScL_disp[n].TM[0] = t_matrix_coeff[n];
    }
    
    cout << "***Building scatterers: done" << endl << endl;
    
    


    // *** Computations ***

    // Construct time-harmonic response
    complex val, k_out;
    Pvec src_wv(1., PI/2., 0., Pvec::Spherical);
    PlaneWave* IW = new PlaneWave(1., src_wv);
    
    for( int s = 0; s < N_src; s++ ){
      
      std::vector< std::vector<complex> > trace_harm(N_rcv, std::vector<complex>(N_freq[s]) );
      std::vector< std::vector<complex> > trace_harm_disp(N_rcv, std::vector<complex>(N_freq[s]) );
      
      
      for( int k = 0; k < N_freq[s]; k++ ){
	cout << "Source idx : " << (s + 1) << "/" << N_src << "  ;  " << "Frequency idx : " << k << "/" << N_freq[s]  << endl;
	
	// Update incoming wave
	src_wv = src_xi[s][k];
	IW->wv = src_wv;
	


	// --- Undisplaced field ---
	
	// Solve linear system (single-scattering; equals r.h.s.)
	solver.RHSInit(IW, ScL, u);  
	
	/*	// Transfer espansion of each scatterer expansion to center
	std::vector<complex> u_trans(local_Index.size());
	for( int n = 0; n < NScat; n ++ ){
	  Pvec p = -ScL[n].getLocation();
	  complex wv = src_wv.r;
	  FMM::PointAndShoot trans(p.r, wv, p.theta, p.phi, &local_Index, &Index, false);
	  std::vector<complex> vec = trans.Apply(u[n]);	
	  
	  for( int i = 0; i < local_Index.size(); i++ )
	    u_trans[i] += vec[i];
	}
	
	// Compute response in far-field
	for( int q = 0; q < N_rcv; q++ ){ 
	  trace_harm[q][k] = 0;
	  for( int i = 0; i < local_Index.size(); i++ ){
	    int l = local_Index(i,0);
	    int m = local_Index(i,1);
	    
	    trace_harm[q][k] += C[s][q][k] * u_trans[i] * 1./(src_wv.r*pow(CI,l+1.)) * gsl_sf_harmonic(l, m, rcv_loc[q].theta, rcv_loc[q].phi);
	    // TODO: exp(CI*src_wv.r*R)/R factor unnecessary right?
	    // exp(CI*src_wv.r*R)/R;

	    //cout << trace_harm[q][k] << endl;
	  }
	  }*/


	for( int q = 0; q < N_rcv; q++ ){
	  trace_harm[q][k] = 0;
	  for( int n = 0; n < NScat; n ++ ){
	    double dot = (rcv_loc[q].x*ScL[n].getLocation().x 
			  + rcv_loc[q].y*ScL[n].getLocation().y 
			  + rcv_loc[q].z*ScL[n].getLocation().z) / rcv_loc[q].r;
	    trace_harm[q][k] += C[s][q][k] * u[n][0] * exp(-CI*src_wv.r*dot);
	  }
	}



	
	// --- Displaced field ---
	
	// Solve linear system (single-scattering; equals r.h.s.)
	solver.RHSInit(IW, ScL_disp, u); 
      

	/*	std::vector<complex> u_trans_disp(local_Index.size());
	for( int n = 0; n < NScat; n ++ ){
	  Pvec p = -ScL_disp[n].getLocation();
	  FMM::PointAndShoot trans(p.r, src_wv.r + 0.*CI, p.theta, p.phi, &local_Index, &Index, false);
	  std::vector<complex> vec = trans.Apply(u[n]);	

	  for( int i = 0; i < local_Index.size(); i++ )
	    u_trans_disp[i] += vec[i];
	}
							

	// Compute response in far-field
	for( int q = 0; q < N_rcv; q++ ){ 
	  trace_harm_disp[q][k] = 0;
	  for( int i = 0; i < local_Index.size(); i++ ){
	    int l = local_Index(i,0);
	    int m = local_Index(i,1);

	    trace_harm_disp[q][k] += C[s][q][k] * u_trans[i] * 1./(src_wv.r*pow(CI,l+1.)) * gsl_sf_harmonic(l, m, rcv_loc[q].theta, rcv_loc[q].phi);
	      // TODO: exp(CI*src_wv.r*R)/R factor unnecessary right?
	      // exp(CI*src_wv.r*R)/R;
	  }
	}*/


							
	// Transfer each scatterer expansion to center
	for( int q = 0; q < N_rcv; q++ ){ 
	  trace_harm_disp[q][k] = 0;
	  for( int n = 0; n < NScat; n ++ ){
	    double dot = (rcv_loc[q].x*ScL_disp[n].getLocation().x 
			  + rcv_loc[q].y*ScL_disp[n].getLocation().y 
			  + rcv_loc[q].z*ScL_disp[n].getLocation().z) / rcv_loc[q].r;
	    
	    trace_harm_disp[q][k] += C[s][q][k] * u[n][0] * exp(-CI*src_wv.r*dot);
	  }
	}
	
							}
      
      
      // TODO: Optimize with fft?
      // Construct time response for each receiver
      for( int q = 0; q < N_rcv; q++ ){ 
	cout << "time response : " << q << endl;
	for( int k = 0; k < N_freq[s]; k++ ){
	  Pvec p = src_xi[s][k];
	  
	  for( int l = 0; l < N_t; l++ ){
	    trace[s][q][l] += trace_harm[q][k] * exp(-CI*p.r*time[q][l]);
	    trace_disp[s][q][l] += trace_harm_disp[q][k] * exp(-CI*p.r*time[q][l]);
	  }
	  
	}	
      }
      
      }					       


    // *** Save data ***
    
    // Recorded signal
    std::string file;
    for( int s = 0; s < N_src; s++){

      std::ostringstream stm ;
      stm << s;
      std::string sgn_filename_ = sgn_filename + "_" + stm.str() + ".csv";
      std::string sgn_filename_d_ = sgn_filename_d + "_" + stm.str() + ".csv";

      ofstream sgn_file(sgn_filename_.c_str(), ios::out);
      ofstream sgn_file_d(sgn_filename_d_.c_str(), ios::out);
      for( int i = 0; i < N_t; i++ ){
	for( int q = 0; q < N_rcv; q++ ){
	  sgn_file << std::setprecision(15) << time[q][i] << " , " << std::real(trace[s][q][i]) << "," << std::imag(trace[s][q][i]) << ",";
	  sgn_file_d << std::setprecision(15) << time[q][i] << " , " << std::real(trace_disp[s][q][i]) << "," << std::imag(trace_disp[s][q][i]) << ",";
	}
	sgn_file << endl;
	sgn_file_d << endl;
      }
      sgn_file.close();
      sgn_file_d.close();
    }


  return 0;
}






#endif
