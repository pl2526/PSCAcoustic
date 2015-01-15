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

#include "./Aux/Eigen/Dense"

inline double T(int n, double a, double b, double r){
  double C = sqrt(2./b);
  if( n == 0 ) C = 1/b;
  return C * Cheb(n, 2./(b-a) * r - (b+a)/(b-a) );
}

inline complex phi(complex k, Pvec& x, Pvec& y){
  Pvec z = x - y;
  complex val = exp(k*CI*z.r)/(4*PI*z.r);
  if( z.r < 1e-10 ) val = 0;
  return val;
}


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
  std::string filename("../AcousticOutputFiles/FFH/HomogenizationCheck_single_A_");
  std::string response_filename("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_single_A_");

  std::string filename_full("../AcousticOutputFiles/FFH/HomogenizationCheck_full_A_");
  std::string response_filename_full("../AcousticOutputFiles/FFH/HomogenizationCheck_farfield_full_A_");



  // *** Full solution (Lippmann-Schwinger) parameters ***
  double C_IN = 0.7;
  double R = 1.;       // "Inside" domain
  int L = 4;             // Max degree expansion for translation at orig
  int K = 9;           // Number of Chebyshev points per cell
  Indexing Index_FF(L);

  // Should be overconstarined
  //assert((L1+1)*(L1+1) >= K*(L2+1)*(L2+1));

  

  // *** Source parameters ***
  std::string src_type("PlaneWave");
  Pvec src_wv(0., 0., std::abs(K_OUT), Pvec::Cartesian);   // Plane wave, eave vector


  // *** Sampling points ***

  int N_S_MAX = 30;
  int N_S = 0;
  std::vector<Pvec> samples;

  for( int i = 0; i < N_S_MAX; i++ )
    samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, 0, Pvec::Spherical) );
  
  for( int i = 0; i < N_S_MAX; i++ )
   samples.push_back( Pvec(1., i/(N_S_MAX-1.)*PI, PI/2., Pvec::Spherical) );
  
  N_S = samples.size();
 



    // *** Initialization of source ***
    PlaneWave* IW = new PlaneWave(1., src_wv);
    






    // Compute quadratures
    std::vector<double> cheb_nodes(2*K);
    std::vector<double> cheb_weights(2*K);
    std::vector<double> gauss_nodes(2*L+1);
    std::vector<double> gauss_weights(2*L+1);
    std::vector<double> trapz_nodes(2*L+1);
    std::vector<double> trapz_weights(2*L+1);
  
    getGaussChebyshevQuad(2*K, cheb_nodes, cheb_weights);
    getGaussLegendreQuad(2*L+1, gauss_nodes, gauss_weights);
    for( int k = 0; k < (2*L+1); k++ ){
      trapz_nodes[k] = k * 2*PI/(2*L+1);
      trapz_weights[k] = 2*PI/(2*L+1);
    }

    std::vector<Pvec> nodes;
    std::vector<double> weights;
    for( int i = 0; i < (int) cheb_nodes.size(); i++ ){
      for( int j = 0; j < (int) gauss_nodes.size(); j++ ){
	for( int k = 0; k < (int) trapz_nodes.size(); k++ ){
	  // R*(cheb_nodes[i]+1.)/2
	  Pvec p( R*(cheb_nodes[i]+1.)/2, PI*(gauss_nodes[j]+1.)/2., trapz_nodes[k], Pvec::Spherical);
	  double w = R/2.*cheb_weights[i] * PI/2.*gauss_weights[j] *trapz_weights[k];

	  nodes.push_back(p);
	  weights.push_back(w);
	}
      }
    }


    // Checking quadrature
    
    complex val = 0;
    //T(0,0,R,nodes[k].r)
    for( int k = 0; k < (int) nodes.size(); k++ )      
      val += weights[k] * T(5,0,R,nodes[k].r) * gsl_sf_harmonic(3, -3, nodes[k].theta, nodes[k].phi) *
	T(5,0,R,nodes[k].r) * std::conj(gsl_sf_harmonic(3, 3, nodes[k].theta, nodes[k].phi)) 
* sin(nodes[k].theta);// * pow(nodes[k].r, (double) m) * pow(nodes[k].theta, (double) n) * pow(nodes[k].phi, (double) o);
    
    double exact = 0.;
    cout << val << " : " << exact << " : " << std::abs( val -  exact) << endl << endl;
    



    // Compute full solution
    Matrix<complex, Dynamic, 1> x;
    Matrix<complex, Dynamic, Dynamic> A; A.resize( K*(L+1)*(L+1), K*(L+1)*(L+1) );    // Linear system (matrix A in Ax = b)
    Matrix<complex, Dynamic, 1> b; b.resize( K*(L+1)*(L+1) );    // Linear system (vector b in Ax = b)

    
    
    // Compute values of convolution
    std::vector< std::vector<complex> > alpha( nodes.size(), std::vector<complex>(K*(L+1)*(L+1)) );
    for( int w = 0; w < (int) nodes.size(); w++ ){	
      cout << "computing : " << w << " / " << nodes.size()  << endl;  

      int j = 0;
      for( int p = 0; p <= L; p++ ){
	for( int q = -p; q <= p; q++){
	  for( int r = 0; r < K; r++){
	    
	    alpha[w][j] = 0;
	    for( int z = 0; z < (int) nodes.size(); z++ ){
	      alpha[w][j] += weights[z] * 
		phi(K_OUT, nodes[w] ,nodes[z]) *
		T(r,0,R,nodes[z].r)*gsl_sf_harmonic(p, q, nodes[z].theta, nodes[z].phi) * 
		( 4*PI*PI/(C_IN*C_IN) -  4*PI*PI/(C_OUT*C_OUT) ) *
		nodes[z].r*nodes[z].r*sin(nodes[z].theta);
	    }

        
	    j++;
	  }
	}
      }
     
    }

     cout << "Computing entries matrix A..." << endl;
      int i = 0;
      for( int p = 0; p <= L; p++ ){
	for( int q = -p; q <= p ; q++){
	  for( int r = 0; r < K; r++){

	    cout << "row : " << i << " / " << K*(L+1)*(L+1)  << endl;
	    
	    // Compute entries of A
	    for( int j = 0; j < K*(L+1)*(L+1); j++ ){

	      A(i,j) = 0;
	      for( int w = 0; w < (int) nodes.size(); w++ ){
		A(i,j) += weights[w] 
		  * T(r,0,R,nodes[w].r)*std::conj(gsl_sf_harmonic(p, q, nodes[w].theta, nodes[w].phi))
		  * alpha[w][j] 
		  *sin(nodes[w].theta);
	      }
	      
	      if( i == j )
		A(i,j) += 1.;

	    }


	    // Compute entries for r.h.s.
	    b(i) = 0;
	    for( int z = 0; z < (int) nodes.size(); z++ ){
	      b(i) += weights[z] *
	        T(r,0,R,nodes[z].r)*std::conj(gsl_sf_harmonic(p, q, nodes[z].theta, nodes[z].phi)) * 
		IW->Evaluate(nodes[z]) *
		(4*PI*PI/(C_OUT*C_OUT) - 4*PI*PI/(C_IN*C_IN)) 
		*sin(nodes[z].theta);
	    }
	    
	    i++;
	  }
	}
      }

  cout << "done" << endl;



  // Solve linear system
  cout << "Solving linear system...";
  cout << endl << "rank : " << A.fullPivHouseholderQr().rank() << " : " << K*(L+1)*(L+1) <<  endl;
  x = A.fullPivHouseholderQr().solve(b);
  cout << "done" << endl;

  
  /*Matrix<complex, Dynamic, 1> z;
  z = A*x;
  for( int i = 0; i < b.size(); i++ )
  cout << "check : " << z(i) << " : " << b(i) << endl;*/
  

  // for( int i = 0; i < x.size(); i++ )
  // cout << x(i) << endl;
  //cout << endl << endl; 

      
      // Compute response
  cout << "Compute response...";
  std::vector<complex> FF( (L+1)*(L+1) );
  i = 0;
  for( int l = 0; l <= L; l++ ){
    for( int m = -l; m <= l; m++){
      FF[i] = 0.;

      int j = 0;
      for( int p = 0; p <= L; p++ ){
	for( int q = -p; q <= p ; q++){
	  for( int r = 0; r < K; r++){
	    if( p == l && m == q )
	      FF[i] +=  x[j] * T(r,0,R,R) / gsl_sf_hankel_1(l, std::abs(K_OUT) * R);
	    j++;
	  }
	}
      }
      
      i++;
    }
  }
      
    
  std::vector<complex> response_LM(samples.size());
  for( int s = 0; s < (int) samples.size(); s++ ){
    for( int i = 0; i < ((L+1)*(L+1)); i++ ){
      int l = Index_FF(i,0);
      int m = Index_FF(i,1);
      
      //cout << FF[i] <<  endl;
      response_LM[s] += 1./(K_OUT*pow(CI, l+0.)) * FF[i] * gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
    }
  }
  cout << "done" << endl;


      
  // Write response
  cout << "   ***Writing response...";
  std::string temp_filename_full =  response_filename_full;
  temp_filename_full.append(".csv");
  char *Filename = (char*) temp_filename_full.c_str();
  ofstream Resfile_full(Filename, ios::out);
  for( int i = 0; i < samples.size(); i++ ){
    Resfile_full << samples[i].theta << ", " << samples[i].phi << ", " 
		 << std::real(response_LM[i]) << ", " << std::imag(response_LM[i]) << endl;
  }
  cout << " done" << endl;
  
      
    
    

  // *** Single sphere ***
  Pvec location(0.,0.,0., Pvec::Cartesian);
  Scatterer scatterer(R, 2*PI/C_IN, K_OUT, 1., location, 0);
  std::vector<complex> vec = IW->coeff;
  scatterer.TM_Apply(vec);



  std::vector<complex> response_S(samples.size());
  for( int s = 0; s < (int) samples.size(); s++ ){
    for( int i = 0; i < ((L+1)*(L+1)); i++ ){
      int l = Index_FF(i,0);
      int m = Index_FF(i,1);

      cout << FF[i] << " : " << vec[i] << endl;

      response_S[s] += 1./(K_OUT*pow(CI, l+0.)) * vec[i] * gsl_sf_harmonic(l, m, samples[s].theta, samples[s].phi);
    }
  }
  cout << "done" << endl;


      // Write response
  cout << "   ***Writing response...";
  std::string temp_filename =  response_filename;
  temp_filename.append(".csv");
  char *Filename2 = (char*) temp_filename.c_str();
  ofstream Resfile(Filename2, ios::out);
  for( int i = 0; i < samples.size(); i++ ){
    Resfile << samples[i].theta << ", " << samples[i].phi << ", " 
		 << std::real(response_S[i]) << ", " << std::imag(response_S[i]) << endl;
  }
  cout << " done" << endl;
  
      
    
    
    
  
  return 0;
    
}




#endif
