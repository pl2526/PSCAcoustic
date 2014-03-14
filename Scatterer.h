#ifndef SCATTERER_H
#define SCATTERER_H

/*
 *  Scatterer.h
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

#include "Constants.h"
#include "Coordinates.h"
#include "Indexing.h"


class Scatterer
{
	
public:

  static Indexing global_index;               // Index for multipole coefficients

  double radius;                              // Radius of the scatterer
  complex  k, k_out;        // Wave vectors
  double rho;                                 // Density of the scatterer
  Pvec location;                              // Location of the scatterer

  // Convention
  // Outer index: degree and order (l,m) of spherical wave functions
  // Inner index: in order from 0 to 9: T_11, T_12, T_13, T_21, T_22, T_23, T_31, T_32, T_33
  std::vector<complex> TM;                   // T-matrix storage

  // Constructors
  Scatterer(){}  // Empty constructor
  Scatterer( double radius_, complex  k_, complex k_out_, double rho_, Pvec& location_ );
  
  // Destructor
  ~Scatterer(){}



  // T-matrix utilities
  static void TMconstruct(std::vector<complex>& TM, complex  k, complex  k_out,  double r, double rho); // Compute entries
  void TM_Apply( std::vector<complex>& incoming);  // Apply T-matrix (in-place)

  // Evaluate VSWF expansion with center c at point p 
  static complex Evaluate(Pvec& p, Scatterer& S, std::vector<complex>& a, SWF_Type type, complex k_out = K_OUT);
  static complex Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type, complex k_out = K_OUT);
  //static complex FF_Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type);
  //static std::vector<complex> Evaluate(Pvec& p, Pvec& c, std::vector<complex>& a, ECC::type Type);
  
  // Accessors
  inline std::vector<complex> getTMatrix() { return TM; }
  inline Pvec getLocation() { return location; }
  inline double getRadius() { return radius; }

  // Copier
  void copy(Scatterer& scat);
  
};



#endif