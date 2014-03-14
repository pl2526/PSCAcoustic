#ifndef IMAGING_CPP
#define IMAGING_CPP

/*
 *  Imaging.cpp
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 11/11/13.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Produce image of wave field
 */


#include "Imaging.h"

void Imaging( std::vector<Scatterer>& ScL, IncWave* IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in){
  std::vector<IncWave*> v(1); v[0] = IW;
  Imaging(ScL, v, filename, u, idx, M, L, u_in);
}

// TODO: receive normal to plane
void Imaging( std::vector<Scatterer>& ScL, std::vector<IncWave*>& IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in)
{
  cout << "   ***Creating images..." << endl;
  // Create a 2D grid onto which to compute field

  std::vector<double> X(2*M+1);
  for( int i = 0; i <= 2*M; i++ )
    X[i] = L * (i-M)/ ((double) M);

  // Compute image and store
  complex v_sc;
  complex v_in;

  // Create files
  std::stringstream idx_str;
  idx_str << idx;

  std::string str_x(filename);
  str_x.append("_xy_");
  str_x.append(idx_str.str());
  str_x.append(".csv");
  char *Str_x = (char*) str_x.c_str();
  ofstream im_x(Str_x, ios::out);

  std::string str_grid = "grid_";
  str_grid.append(idx_str.str());
  str_grid.append(".csv");
  char *Str_grid = (char*) str_grid.c_str();
  ofstream grid(Str_grid, ios::out);



  // TODO: Transfer expansion from all transducers to the origin to obtain single expansion
  for( int i = 0; i < X.size(); i++ ){
    //cout << "row : " << i << endl;
    for( int j = 0; j < X.size(); j++ ){
      //cout << "col : " << j << endl;
      Pvec p(0.,X[i],X[j],Pvec::Cartesian);

      bool in = false;
      int in_idx = 0;
      Pvec q,s;
      for( int n = 0; n < NScat; n++){
	q = ScL[n].getLocation();
	s = p - q;
	if( s.r <= RADIUS ){
	  in = true;
	  in_idx = n;
	  break;
	}
      }

      // Record location of gridpoints
      grid << X[i] << "," << X[j] << "," << 0. << endl;
      
      // Evaluate scattered field at point

      v_sc = 0.;
      v_in = 0.;
      if( in && ~(u_in == NULL) ){
	v_in = Scatterer::Evaluate(p, ScL[in_idx], (*u_in)[in_idx], Bessel);

      } else { // If image point is oustide all scatterers
	// Scattered wave
	v_sc = Scatterer::Evaluate(p, ScL, u, Hankel);
	
	// Incoming wave
	for( int k = 0; k < (int) IW.size(); k++)
	  v_in = IW[k]->Evaluate(p);
      }
      
      im_x << std::setprecision(15) << std::real(v_sc + v_in) << "," << std::setprecision(15) << std::imag(v_sc + v_in) << endl;
    }
  }


  // Close files
  im_x.close();
  grid.close();
  cout << "   ***Creating images: done" << endl;  

  return;
}







#endif
