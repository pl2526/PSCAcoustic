/*
 *  Scatterer.h
 *  PSCAcoustic
 *
 *  Objects describing scatterers.
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
