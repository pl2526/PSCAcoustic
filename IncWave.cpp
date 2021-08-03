/*
 *  IncWave.cpp
 *  PSCAcoustic
 *
 *  Objects for incoming wave field.
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

#ifndef INCWAVE_CPP
#define INCWAVE_CPP

#include "IncWave.h"

Indexing IncWave::global_index(LMAX);

// Spherical wave

// Constructor
SphericalWave::SphericalWave(complex k_out_, complex strength_, const Pvec& loc_) : IncWave(), strength(strength_), loc(loc_)
    {
      k_out = k_out_;
      SW_index = new Indexing(0);
      coeff.resize(1);
      coeff[0] = strength*k_out*CI/sqrt(4*PI);
    }

void SphericalWave::Translate(Pvec& c, std::vector<complex>& T_c)
{ 
  assert( (int) T_c.size() == (int) global_index.size() );
  Pvec vec = c-loc;
  Pvec w = -loc;
  
  // Point source
  if( vec.r > 1e-11  && vec.r < 1e6  ){
    FMM::PointAndShoot SW_Transfer(vec.r, real(k_out), vec.theta, vec.phi, &global_index, SW_index, false);
    // TODO : TRANSFER SHOULD BE DONE IN-PLACE
    std::vector<complex> temp = coeff;
    T_c = SW_Transfer.Apply(temp, FMM::FORWARD);
  } else if( vec.r <= 1e-11 ) {
    T_c.resize(global_index.size());
    T_c[0] = coeff[0];
    for( int i = 1; i < global_index.size(); i++ )
      T_c[i] = 0;
  } else if( vec.r >= 1e6 ) {
    for( int i = 0; i < global_index.size(); i++ )
      T_c[i] = 4*PI *strength * exp(k_out*CI*loc.r)/(4*PI) * 
	std::conj( gsl_sf_harmonic(global_index(i,0), global_index(i,1), w.theta, w.phi) ) * 
	exp(-k_out*CI*(c.x*(loc.x/loc.r) + c.y*(loc.y/loc.r) + c.z*(loc.z/loc.r)) ) * pow(CI,(double) global_index(i,0));
  }
  
  return;
}


complex SphericalWave::Evaluate(Pvec& p){
  complex val = 0;
  Pvec v = p-loc;

  if( v.r < 1e6 ){
    val = strength*1./(4*PI*v.r)*exp(k_out*CI*v.r);
    //cout << "val : " << K_OUT/(4*PI)*CI * gsl_sf_hankel_1(0, std::real(k_out)*v.r) << " : " << 1./(4*PI*v.r)*exp(k_out*CI*v.r) << " : " <<  std::abs( K_OUT/(4*PI)*CI * gsl_sf_hankel_1(0, std::real(k_out)*v.r) - 1./(4*PI*v.r)*exp(k_out*CI*v.r) ) << endl;
  }
  else if( v.r >= 1e6 ){
    double r;
    ( loc.r > p.r ) ? r = loc.r : r = p.r;
   
    val = strength * exp(k_out*CI*r)/(4*PI) 
      * exp(-k_out*CI*(p.x*(loc.x/r) + p.y*(loc.y/r) + p.z*(loc.z/r)) );

  }

  return val;
}




// Plane wave
// The wave vector wv should NOT be normalized
PlaneWave::PlaneWave(complex amp_, Pvec wv_) : IncWave(), amp(amp_), wv(wv_)
{
  k_out = wv.r;
  coeff.resize(global_index.size());
  for ( int i = 0 ; i < global_index.size() ; i++ )
    coeff[i] = 4*PI * amp * pow(CI,(double) global_index(i,0)) * std::conj( gsl_sf_harmonic(global_index(i,0), global_index(i,1), wv.theta, wv.phi) );
  
}

void PlaneWave::Translate( Pvec& c, std::vector<complex>& T_c)
{ 
  assert( (int) T_c.size() == (int) global_index.size() );

  for( int k = 0; k < global_index.size(); k++)
    T_c[k] = coeff[k] * exp( (c.x*wv.x + c.y*wv.y + c.z*wv.z ) * CI );
  
  return;
}

complex PlaneWave::Evaluate(Pvec& p){
  return  amp * exp( CI * (p.x*wv.x + p.y*wv.y + p.z*wv.z ) );
}



#endif
