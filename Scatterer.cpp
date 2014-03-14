#ifndef SCATTERER_CPP
#define SCATTERER_CPP
/*
 *  Scatterer.cpp
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  This file deals with everything pertaining to a given scatterer:
 *  among other things, its T-matrix and local multipole expansion.
 *
 *  A scatterer object has physical properties, a location in space, a T-matrix
 *  that allows to compute its scattered field given the incostd::ming fiels and a 
 *  scattered expansion.
 *
 */

#include <complex>  
#include <vector>  
#include <assert.h>

typedef std::complex<double> complex;

#include "Scatterer.h"
#include "SpecialFunctions.h"



// TODO: should not be creating different versions of Index...
Indexing Scatterer::global_index(LMAX);

// Constructor
Scatterer::Scatterer( double radius_, complex  k_, complex k_out_, double rho_, Pvec& location_ ) : radius(radius_), k(k_), k_out(k_out_), rho(rho_), location(location_)
{
  TM.resize(LMAX+1);
  TMconstruct(TM, k, k_out, radius, rho);  // Construct T-matrix
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
static complex  Mie(complex  k, complex  k_out, double radius, double rho, int l, int m = 0){
    complex  mie;

    // TODO : Check rh_out
    complex  arg_in_C = k * radius;
    complex  arg_out_C = k_out * radius;
    complex  gamma = k/k_out * RHO_OUT/rho;
    
    // Penetrable sphere (w/ damping)
    complex  numerator = gamma * Amos_sf_bessel_jl(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C)
      - Amos_sf_bessel_jl_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C);
    complex  denominator = Amos_sf_hankel_l_prime(l, arg_out_C) * Amos_sf_bessel_jl(l, arg_in_C)
      - gamma * Amos_sf_hankel_1(l, arg_out_C) * Amos_sf_bessel_jl_prime(l, arg_in_C);
    
    mie = numerator / denominator;
  

    //Sound-hard sphere (Martin, p.132)
    //mie = -Amos_sf_bessel_jl_prime(l, arg_out_C) / Amos_sf_hankel_l_prime(l, arg_out_C);

    //Sound-soft sphere
    //mie = -Amos_sf_bessel_jl(l, arg_out_C) / Amos_sf_hankel_1(l, arg_out_C);


    //Static
    /*  if( l == 0 )
      mie = -Amos_sf_bessel_jl_prime(l, arg_out_C) / Amos_sf_hankel_l_prime(l, arg_out_C);
    else
      mie = -Amos_sf_bessel_jl(l, arg_out_C) / Amos_sf_hankel_1(l, arg_out_C);
    */

    return mie;
    
  }

void Scatterer::TMconstruct(std::vector<complex>& TM, complex  k, complex  k_out,  double r, double rho){
  
  for ( int l=0 ; l < (int) TM.size() ; l++)
    {
      if( l != 0 && std::abs(TM[l-1]) < 1e-20 ){
	TM[l] = 0.; 
      }else{
	TM[l] = Mie(k, k_out, r, rho, l);
      }
      //cout << "T-matrix : " << TM[l] << endl;
    }

  return;
}



complex Scatterer::Evaluate(Pvec& p, Scatterer& S, std::vector<complex>& a, SWF_Type type, complex k_out){
  
  complex val = 0;
  Pvec q = p-S.location;
  Pvec loc = S.location;
  
  if( type == Bessel ){
    
    for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ ){
      int l = global_index(i,0);
      int m = global_index(i,1);
      
      val += a[i] * gsl_sf_harmonic(l, m, q.theta, q.phi) * Amos_sf_bessel_jl(l, K*q.r);
    }
  }else if( type == Hankel ) {
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

complex Scatterer::Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type, complex k_out){
  assert( ScL.size() == a.size() );
  
  complex val = 0;
  for( int n = 0; n < ScL.size(); n++ )
    val += Scatterer::Evaluate(p, ScL[n], a[n], type, k_out);

  return val;
}

/*complex Scatterer::FF_Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type){
  assert( ScL.size() == a.size() );
  
  complex val = 0;
  for( int n = 0; n < ScL.size(); n++ ){
    Pvec loc = ScL[n].location;
    Pvec q = p-loc;

      for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ ){
	int l = global_index(i,0);
	int m = global_index(i,1);
	
	val += 1./(K_OUT*pow(CI, l+1.)) * a[n][i] * gsl_sf_harmonic(l, m, q.theta, q.phi);
      }
      
      val *= exp(-k_out*CI*(p.x*(loc.x/loc.r) + p.y*(loc.y/loc.r) + p.z*(loc.z/loc.r)) )
	* exp(k_out*CI*loc.r)/loc.r;
  }

  return val;
}
*/



#endif

