/*
 *  Imaging.cpp
 *  PSCAcoustic
 *
 *  Produces image of wave field.
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


#ifndef IMAGING_CPP
#define IMAGING_CPP

#include "Imaging.h"

void Imaging( std::vector<Scatterer>& ScL, IncWave* IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in, complex k_out, complex k_in){
  std::vector<IncWave*> v(1); v[0] = IW;
  Imaging(ScL, v, filename, u, idx, M, L, u_in, k_out, k_in);
}

// TODO: receive normal to plane
void Imaging( std::vector<Scatterer>& ScL, std::vector<IncWave*>& IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in, complex k_out, complex k_in)
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

  std::string str_Re(filename + "Re_" + idx_str.str() + ".csv");
  std::string str_Im(filename + "Im_" + idx_str.str() + ".csv");

  char *Str_Re = (char*) str_Re.c_str();
  ofstream im_Re(Str_Re, ios::out);
  char *Str_Im = (char*) str_Im.c_str();
  ofstream im_Im(Str_Im, ios::out);


  // TODO: Transfer expansion from all transducers to the origin to obtain single expansion
  for( int i = 0; i < X.size(); i++ ){
    //cout << "row : " << i << " / " << X.size() << endl;
    for( int j = 0; j < X.size(); j++ ){
      //cout << "col : " << j << " / " << X.size() << endl;
      Pvec p(X[i],X[j],0.,Pvec::Cartesian);

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

      // Evaluate scattered field at point
      v_sc = 0.;
      v_in = 0.;
      if( in ){

	v_in = Scatterer::Evaluate(p, ScL[in_idx], (*u_in)[in_idx], Bessel, k_in);

      } else { // If image point is oustide all scatterers
	// Scattered wave
	v_sc = Scatterer::Evaluate(p, ScL, u, Hankel, k_out);
	
	// Incoming wave
	for( int k = 0; k < (int) IW.size(); k++)
	  v_in += IW[k]->Evaluate(p);
      }
      
      im_Re << std::setprecision(15) << std::real(v_sc + v_in) << "," ;
      im_Im << std::setprecision(15) << std::imag(v_sc + v_in) << "," ;
    }

    im_Re << endl;
    im_Im << endl;
  }


  // Close files
  im_Re.close();
  im_Im.close();
  cout << "   ***Creating images: done" << endl;  

  return;
}







#endif
