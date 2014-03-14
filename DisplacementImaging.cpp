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

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "Imaging.h"
#include "PSCEnv.h"
#include "Solver.h"
#include "Finalize.h"

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
  double x_disp = amp * 0.;//exp(-pow(eps, delta)/2. * p.x*p.x);
  double y_disp = 1e-5;//amp * 1.;
  double z_disp = amp * 0;
  Pvec disp(x_disp, y_disp, z_disp, Pvec::Cartesian);
  return disp;
}


// Windows for curvelets
double W_r(double scale, double r){
    double center_r = 5./(4.*scale);
    double sigma_r = -pow(3./(4.*scale),2.)/log(1e-5);
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
  std::string src_type("PW");
  double mu = 1;               // Mean of scattering coefficients
  double a = 1.;                // amplitude of scattering coefficients
  double src_scale = 1e-2;      // Scale of the curvelet
  double delta = pow(src_scale, 1.5);
  double disp_amp_scale =  1e-6;//pow(src_scale, 1.+delta);           // Scale of the amplitude of the displacement




  // *** File names ***
  std::string filename("DI_SingleScattering_mu=0_HH");
  std::string freq_filename("DI_SingleScattering_mu=0_freq_HH");
  std::string extra_filename("DI_SingleScattering_mu=0_extra_HH");
  std::string disp_filename("DI_SingleScattering_mu=0_HH_disp");
  std::string rcv_filename("DI_SingleScattering_mu=0_rcv_HH");
  std::string image_filename("DI_SingleScattering_mu=0_image_HH");

  std::string filename_d(filename); filename_d.append("_disp");
  std::string freq_filename_d(freq_filename); freq_filename_d.append("_disp.csv");
  std::string extra_filename_d(extra_filename); extra_filename_d.append("_disp.csv");
  std::string disp_filename_d(disp_filename); disp_filename_d.append("_disp.csv");
  std::string rcv_filename_d(rcv_filename); rcv_filename_d.append("_disp.csv");
  std::string image_filename_d(image_filename); image_filename_d.append("_disp.csv");

  freq_filename.append(".csv");
  extra_filename.append(".csv");
  disp_filename.append(".csv");
  rcv_filename.append(".csv");
  image_filename.append(".csv");






  // *** Image parameters ***
  bool produce_image = false;
  int M = 35;     // Number of pixels: 2*M+1 from -M to M
  double L;   // size of image [-L,L]x[-L,L]




  // *** Sources' parameters ***
  Pvec src_loc(0.,0.,0.,Pvec::Cartesian);
  Pvec src_dir(1.,PI/2.,PI/2,Pvec::Spherical);    // Curvelet direction
 



// *** Receivers' locations ***
  // TODO: currently only works for 1 receiver
  std::vector<Pvec> rcv_loc(1);
  rcv_loc[0] =  Pvec(1.,0.1,0.,Pvec::Cartesian);
  double max_rcv_loc_r = 0;
  for( int i = 0; i < rcv_loc.size(); i++ )
    max_rcv_loc_r = std::max(max_rcv_loc_r, rcv_loc[i].r);





  // *** Displacement parameters *** 
  ofstream disp_file((char*) disp_filename.c_str(), ios::out);

  int N_d_x = 200;
  int N_d_y = 100;
  double L_d_x = 2;
  double L_d_y = 1;
  for( int i = 0; i <= N_d_x; i++ ){
    for( int j = -N_d_y; j <= N_d_y; j++ ){
      Pvec loc(i*L_d_x/N_d_x, j*L_d_y/N_d_y, 0., Pvec::Cartesian);
      Pvec disp = Displacement(loc, disp_amp_scale, src_scale, delta);
      disp_file << loc.x << "," << loc.y << "," << disp.x <<  "," << disp.y << "," << disp.z << "," << endl;
    }
    disp_file << endl;
  }




  // *** Frequency content of curvelet ***
  double R = 0.15;                                 // Radius of spherical culster within which scatterers lie

  // xi_r corresponds to radius
  // xi_p corresponds to phi
  double xi_r_min = 1./(2.*src_scale);
  double xi_r_max = 2./src_scale;
  double xi_p_min = src_dir.phi - 1.1*sqrt(src_scale);
  double xi_p_max = src_dir.phi + 1.1*sqrt(src_scale);

  // Period considered for imaging purposes
  double T = 1.2*(max_rcv_loc_r + 2*R);          // Period
  int N_t = 15*T / (4*src_scale);
  cout << "T : " << T << ", N_t : " << N_t << endl;
  double src_rate = 1./(2.*T);  
  //double src_rate_r = 1./(2.*T);                 // Sampling rate in radial direction
  //double src_rate_p = src_rate_r / xi_r_max;     // Sampling rate in angular direction
  std::vector<double> time(N_t);
  for( int i = 0; i < N_t; i++ )
    time[i] = ((double) i)/N_t * T;
 

  /*  // Sampling in radial direction
  int src_N_r = std::ceil((xi_r_max - xi_r_min) / src_rate_r);
  std::vector<double> src_xi_r(src_N_r+1);
  for(int i = 0; i <= src_N_r; i++ )
    src_xi_r[i] = xi_r_min + ((double) i)/src_N_r * (xi_r_max - xi_r_min);


  // Sampling in angular direction
  int src_N_p = std::ceil((xi_p_max - xi_p_min) / src_rate_p);
  std::vector<double> src_xi_p(src_N_p+1);
 for(int i = 0; i <= src_N_p; i++ )
 src_xi_p[i] = xi_p_min + ((double) i)/src_N_p * (xi_p_max - xi_p_min);*/

  // Sampling in x direction
  int src_N_y = std::ceil((xi_r_max - xi_r_min) / src_rate);
  std::vector<double> src_xi_y(src_N_y+1);
  for(int i = 0; i <= src_N_y; i++ ){
    src_xi_y[i] = xi_r_min + ((double) i)/src_N_y * (xi_r_max - xi_r_min);
    //cout << " y : " <<  src_xi_y[i] << endl;
  }


  // Sampling in y direction
  double x_min = -xi_r_max*sin(xi_p_max);// - src_dir.phi);
double x_max = xi_r_max*sin(xi_p_max);// - src_dir.phi);
  int src_N_x = 2*std::ceil( (x_max - x_min) / (2*src_rate));
  std::vector<double> src_xi_x(src_N_x+1);
  for(int i = 0; i <= src_N_x; i++ ){
    src_xi_x[i] = x_min + ((double) i)/src_N_x * (x_max - x_min);
    cout << " x : " <<  src_xi_x[i] << endl;
  }

  //cout << xi_r_max << " : " <<sin(xi_p_max) << " : " << sin(xi_p_min) << " : " << src_rate << " : " << endl;	        
 cout << "Number of samples : x : " << src_N_x << " : y : " << src_N_y << endl;



 /*
 // TODO: Clean here
  std::vector< std::vector<double> > src_amp(src_N_r+1, std::vector<double>(src_N_p+1));
  ofstream freq_file((char*) freq_filename.c_str(), ios::out);
  for(int i = 0; i < src_N_r; i++ ){
    for(int j = 0; j < src_N_p; j++ ){
      //freq_file << src_xi_r[i] << "," << src_xi_p[j] << endl;
      src_amp[i][j] = 1.;
    }
  }
  freq_file.close();
 */



  // *** Save extra information to file ***
  // Save extra info to file (info that's not used in the standard cases)
  char *Filename = (char*) extra_filename.c_str();
  ofstream extra_file(Filename, ios::out);
  extra_file << "Number of freq samples in x: ," << src_N_x << endl;
  extra_file << "Number of freq samples in y: ," << src_N_y << endl;
  extra_file << "Freq sampling rate (radius): ," << src_rate << endl;
  extra_file << "Freq sampling rate (phi): ," << src_rate << endl;
  extra_file << "Curvelet scale: ," << src_scale << endl;
  extra_file << "Curvelet direction: ," << src_dir.x << " , "<<  src_dir.y << " , " << src_dir.z << " , " << endl;
  extra_file << "Image window center: ," << N_d_x << endl;
  extra_file << "Image window N_x pixels: ," << N_d_x << endl;
  extra_file << "Image window N_y pixels: ," << N_d_y << endl;
  extra_file << "Image window L_x size: ," << L_d_x << endl;
  extra_file << "Image window L_y size: ," << L_d_y << endl;








  // ---- Computations -----

  // *** Construct scatterers ***
  cout << "***Building scatterers..." << endl;
  std::vector<Scatterer> ScL;
  std::vector<Scatterer> ScL_disp;
  double d_min = 0.;
  Pvec center(0.,1.1*R,0.,Pvec::Cartesian);
  //std::vector<Pvec> scat_loc(1); scat_loc[0] = Pvec(0.,0.,0.,Pvec::Cartesian);
  std::vector<Pvec> scat_loc = RandSphericalXY(R, RADIUS, d_min, center, NScat);
  assert(scat_loc.size() == NScat);
  std::vector<complex> t_matrix_coeff(NScat); // Nonscales  
  for( int i = 0; i < NScat; i++ ){
    
    Pvec loc(scat_loc[i].x, scat_loc[i].y, 0., Pvec::Cartesian);
    Pvec loc_disp = loc + Displacement(loc, disp_amp_scale, src_scale, delta);
    t_matrix_coeff[i] = 1.;//(mu + 0.1*mu*(2.*drand48()-1.)); // The Mie coefficients are overriden
    
    Scatterer scatterer( RADIUS, K, K_OUT, RHO, loc);
    //scatterer.TM[0] = t_matrix_coeff[i];
    ScL.push_back(scatterer);

    Scatterer scatterer_disp( RADIUS, K, K_OUT, RHO, loc_disp);
    //scatterer_disp.TM[0] = t_matrix_coeff[i];
    ScL_disp.push_back(scatterer_disp);

  }
  cout << "***Building scatterers: done" << endl << endl;
  


  // *** Compute solutions for each source and receiver locations ***

  std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
  double res, rel_res, cond;
  Proc_Idx = 0;
  std::vector< std::vector<complex> > trace(rcv_loc.size(), std::vector<complex>(N_t));
  std::vector< std::vector<complex> > trace_disp(rcv_loc.size(), std::vector<complex>(N_t));
  //std::vector< std::vector<complex> > trace_harm(src_N_r, std::vector<complex>(src_N_p));
  //std::vector< std::vector<complex> > trace_harm_disp(src_N_r, std::vector<complex>(src_N_p));
  std::vector< std::vector<complex> > trace_harm(src_N_x, std::vector<complex>(src_N_y));
  std::vector< std::vector<complex> > trace_harm_disp(src_N_x, std::vector<complex>(src_N_y));
  for( int i = 0; i < rcv_loc.size(); i++ ){
    for( int j = 0; j < N_t; j++ ){
      trace[i][j] = 0.;
      trace_disp[i][j] = 0.;
    }
  }



  // Construct time-harmonic response
  complex val, k_out, T_coeff;
  Pvec src_wv(1., PI/2., 0., Pvec::Spherical);
  PlaneWave* IW = new PlaneWave(1., src_wv);
  for( int i = 0; i < src_N_x; i++ ){
    cout << i << endl;
    

    for( int j = 0; j < src_N_y; j++ ){
      //Proc_Idx = i*src_N_p + j;
      //double C = W_r(src_scale, src_xi_r[i]) * W_p(src_scale, src_dir.phi, src_xi_p[j]) / (src_N_r * src_N_p);


      // Initialization of source (r.h.s.)
      //Pvec src_wv(src_xi_r[i], PI/2., src_xi_p[j], Pvec::Spherical);
      Pvec src_wv(src_xi_x[i], src_xi_y[j], 0., Pvec::Cartesian);
      IW->wv = src_wv;
      //cout << src_xi_x[i] << " : " << src_xi_y[i] << endl;
      //cout <<  src_wv.phi << " : " << src_dir.phi << " : " <<  src_wv.phi - src_dir.phi << endl;
      double C = W_r(src_scale, src_wv.r) * W_p(src_scale, src_dir.phi, src_wv.phi) / (src_N_x*src_N_y);
      //cout << "C : " << C << endl;
      for( int n = 0; n < NScat; n++ ){
	T_coeff = src_wv.r*src_wv.r * t_matrix_coeff[n]; // The Mie coefficients are overriden
	ScL[n].TM[0] = T_coeff;
	ScL_disp[n].TM[0] = T_coeff;
      }
	

	// --- Undisplaced field ---
	
	// Solve linear system (single-scattering; equals r.h.s.)
	//cout << "***Initializing right-hand side..." << endl;
	solver.RHSInit(IW, ScL, u);  
	//cout << "***Initializing right-hand side: done" << endl;
	
	// Write information about problem to file
	/* cout << "   ***Writing to file..." << endl;
	   Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename);
	   Write_Source(Proc_Idx, IW, IncWave::Pt, filename);
	   Write_Location(Proc_Idx, ScL, filename);
	   Write_Solution( Proc_Idx, Index, u, filename);
	   cout << "   ***Writing to file: done" << endl << endl;
	*/

	// Record signal at receivers location
	trace_harm[i][j] = C *IW->Evaluate(rcv_loc[0]);// (Scatterer::Evaluate(rcv_loc[0], ScL, u, Hankel, src_wv.r));// + IW->Evaluate(rcv_loc[0]));//
	  
       
	

	
	// Produce image
	//if( produce_image )
	//  Imaging(ScL, IW, image_filename, u, Proc_Idx, M, L);



	// --- Displaced field ---
	
	// Solve linear system (single-scattering; equals r.h.s.)
	//cout << "***Initializing right-hand side..." << endl;
	solver.RHSInit(IW, ScL_disp, u);  
	//cout << "***Initializing right-hand side: done" << endl;
	
	// Write information about problem to file
	/* cout << "   ***Writing to file..." << endl;
	   Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename);
	   Write_Source(Proc_Idx, IW, IncWave::Pt, filename);
	   Write_Location(Proc_Idx, ScL, filename);
	   Write_Solution( Proc_Idx, Index, u, filename);
	   cout << "   ***Writing to file: done" << endl << endl;
	*/

	// Record signal at receivers location
	trace_harm_disp[i][j] = IW->Evaluate(rcv_loc[0]);// C * (Scatterer::Evaluate(rcv_loc[0], ScL_disp, u, Hankel, src_wv.r));// + IW->Evaluate(rcv_loc[0]));
	

    }
  }


  // Construct time response
  cout << "time response" << endl;
  for( int i = 0; i < src_N_x; i++ ){
    cout << i << endl;
    for( int j = 0; j < src_N_y; j++ ){
      Pvec p(src_xi_x[i], src_xi_y[j], 0., Pvec::Cartesian);
      for( int k = 0; k < rcv_loc.size(); k++ ){
	for( int l = 0; l < N_t; l++ ){
	  //trace[k][l] += trace_harm[i][j] * src_xi_r[i] * exp(-CI*src_xi_r[i]*time[l]);
	  //trace_disp[k][l] += trace_harm_disp[i][j] * src_xi_r[i] * exp(-CI*src_xi_r[i]*time[l]);

	  trace[k][l] += trace_harm[i][j] * exp(-CI*p.r*time[l]);
	  trace_disp[k][l] += trace_harm_disp[i][j] * exp(-CI*p.r*time[l]);
	}
      }
    }
  }
	  



  // Save data
  ofstream rcv_file(rcv_filename.c_str(), ios::out);
  ofstream rcv_file_d(rcv_filename_d.c_str(), ios::out);
  for( int k = 0; k < rcv_loc.size(); k++ ){
    rcv_file << rcv_loc[k].x << "," <<  rcv_loc[k].y << ",";
    rcv_file_d << rcv_loc[k].x << "," <<  rcv_loc[k].y << ",";
  }
  rcv_file << endl;
  rcv_file_d << endl;
 

  
  for( int i = 0; i < N_t; i++ ){
    for( int k = 0; k < rcv_loc.size(); k++ ){
      rcv_file << time[i] << " , " << std::real(trace[k][i]) << "," << std::imag(trace[k][i]) << ",";
      rcv_file_d << time[i] << " , " << std::real(trace_disp[k][i]) << "," << std::imag(trace_disp[k][i]) << ",";
    }
    rcv_file << endl;
    rcv_file_d << endl;
  }
  rcv_file.close();
  rcv_file_d.close();



  return 0;
}






#endif
