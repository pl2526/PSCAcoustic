/*
 *  Scatterer.h
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


#ifndef SCATTERER_H
#define SCATTERER_H

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
  Scatterer( double radius_, complex  k_, complex k_out_, double rho_, Pvec& location_, int flag = 0);
  
  // Destructor
  ~Scatterer(){}



  // T-matrix utilities
  static void TMconstruct(std::vector<complex>& TM, complex  k, complex  k_out,  double r, double rho, int flag); // Compute entries
  void TM_Apply( std::vector<complex>& incoming);  // Apply T-matrix (in-place)
  static complex  Mie(complex  k, complex  k_out, double radius, double rho, int flag, int l, int m = 0);
  static complex Mie_in(complex  k, complex  k_out, double radius, double rho, double rho_out, int l); // T-matrix coefficients for inside field

  // Evaluate VSWF expansion with center c at point p 
  static complex Evaluate(Pvec& p, Scatterer& S, std::vector<complex>& a, SWF_Type type, complex wave_number);
  static complex Evaluate(Pvec& p, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& a, SWF_Type type, complex wave_number);
  
  // Accessors
  inline std::vector<complex> getTMatrix() { return TM; }
  inline Pvec getLocation() { return location; }
  inline double getRadius() { return radius; }

  // Copier
  void copy(Scatterer& scat);
  
};



#endif
