/*
 *  Scatterer.cpp
 *  PSCAcoustic
 *
 *  Objects describing scatterers.
 *
 *
 *  Copyright (C) 2014 Pierre-David Letourneau
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef SCATTERER_CPP
#define SCATTERER_CPP

#include <complex>  
#include <vector>  
#include <assert.h>

typedef std::complex<double> complex;

#include "Scatterer.h"
#include "SpecialFunctions.h"



// TODO: should not be creating different versions of Index...
Indexing Scatterer::global_index(LMAX);

// Constructor
Scatterer::Scatterer( double radius_, complex  k_, complex k_out_, double rho_, Pvec& location_, int flag) : radius(radius_), k(k_), k_out(k_out_), rho(rho_), location(location_)
{
  TM.resize(LMAX+1);
  TMconstruct(TM, k, k_out, radius, rho, flag);  // Construct T-matrix
} 


// Application of T-matrix (in-place)
// convention: 0 -> curl; 1 -> curl_curl; 2 -> gradient; 
void Scatterer::TM_Apply( std::vector<complex>& incoming){
  
  // Sanity checks
  assert( (int) incoming.size() == global_index.size() ); 
  std::vector<complex> v(incoming.size());

  // Apply T-matrix in VSWF coordinate system
  for( int i = 0; i < global_index.size(); i++ ){
    int l = global_index(i,0);

    //cout << "TM : " << TM[l] << endl;
    v[i] = TM[l] * incoming[i];
  }

  incoming = v;

  return;
}



// *** Auxiliary functions ***

void Scatterer::copy(Scatterer& scat){
  radius = scat.radius;
  k = scat.k;
  k_out = scat.k_out;
  rho = scat.rho;
  location = scat.location;
  TM = scat.TM;
  
  return;
}



// Compute entries of T-matrix
complex  Scatterer::Mie(complex  k, complex  k_out, double radius, double rho, int flag, int l, int m){
    complex  mie;

    //cout << "flag : " << flag << endl;

    // TODO : Check rh_out
    complex  arg_in_C = k * radius;
    complex  arg_out_C = k_out * radius;
    complex  gamma = k/k_out * RHO_OUT/rho;

    if( flag == 0 ){
      // Penetrable sphere (w/o damping)
      complex  numerator = gamma * Amos_sf_bessel_jl(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C)
	- Amos_sf_bessel_jl_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C);
      complex  denominator = Amos_sf_hankel_l_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C)
	- gamma * Amos_sf_hankel_1(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C);
      
      mie = numerator / denominator;   
    } else if( flag == 1 ){
      //Sound-hard sphere (Martin, p.132)
      mie = -Amos_sf_bessel_jl_prime(l, arg_out_C) / Amos_sf_hankel_l_prime(l, arg_out_C);
    } else if (flag == 2 ){
      //Sound-soft sphere
      mie = -Amos_sf_bessel_jl(l, arg_out_C) / Amos_sf_hankel_1(l, arg_out_C);
    } else if (flag == 3 ){
      //Static
      if( l == 0 )
	mie = -Amos_sf_bessel_jl_prime(l, arg_out_C) / Amos_sf_hankel_l_prime(l, arg_out_C);
      else
	mie = -Amos_sf_bessel_jl(l, arg_out_C) / Amos_sf_hankel_1(l, arg_out_C);
    }

    return mie;
    
  }


complex Scatterer::Mie_in(complex  k, complex  k_out, double radius, double rho, double rho_out, int l){
  double kappa = rho / rho_out;
  complex mie_in = -k_out * radius*radius * CI * ( k_out*Amos_sf_hankel_l_prime(l, k_out*radius) * Amos_sf_bessel_jl(l, k*radius)
						   - k/kappa * Amos_sf_hankel_1(l, k_out*radius) * Amos_sf_bessel_jl_prime(l, k*radius) );
  
  return mie_in;  
}

void Scatterer::TMconstruct(std::vector<complex>& TM, complex  k, complex  k_out,  double r, double rho, int flag){
  
  for ( int l=0 ; l < (int) TM.size() ; l++)
    {
      if( l != 0 && std::abs(TM[l-1]) < 1e-20 ){
	TM[l] = 0.; 
      }else{
	TM[l] = Mie(k, k_out, r, rho, flag, l);
      }
      //cout << "T-matrix : " << TM[l] << endl;
    }

  return;
}



complex Scatterer::Evaluate(Pvec& p, Scatterer& S, std::vector<complex>& a, SWF_Type type, complex wave_number){
  
  complex val = 0;
  Pvec q = p-S.location;
  Pvec loc = S.location;
  
  if( type == Bessel ){   
    complex k = wave_number;

    for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ ){
      int l = global_index(i,0);
      int m = global_index(i,1);
      
      val += a[i] * gsl_sf_harmonic(l, m, q.theta, q.phi) * Amos_sf_bessel_jl(l, k*q.r);
    }
  }else if( type == Hankel ) {
    complex k_out = wave_number;

    for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ ){
      int l = global_index(i,0);
      int m = global_index(i,1);
      
      if( q.r < 1e5 )
	val += a[i] * gsl_sf_harmonic(l, m, q.theta, q.phi) * Amos_sf_hankel_1(l, k_out*q.r);
      else { // Normalized far-field
	double r = std::max(loc.r, p.r);
	
	val += a[i] * gsl_sf_harmonic(l, m, q.theta, q.phi) * 1./(k_out*pow(CI,l+1.)) * 
	  exp(-(p.x*(loc.x/r) + p.y*(loc.y/r) + p.z*(loc.z/r))*k_out*CI) * exp(r*k_out*CI);
      }
      
    } 
    

  }

  return val;
}

complex Scatterer::Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type, complex wave_number){
  assert( ScL.size() == a.size() );
  
  complex val = 0;
  for( int n = 0; n < ScL.size(); n++ )
    val += Scatterer::Evaluate(p, ScL[n], a[n], type, wave_number);

  return val;
}



#endif


