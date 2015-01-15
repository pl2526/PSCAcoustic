/*
 *  IncWave.h
 *  PSCAcoustic
 *
 *  Objects for incoming wave field.
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

#ifndef INCWAVE_H
#define INCWAVE_H

#include "General.h"
#include "./FMMPS/TransferUtils.h"
#include "Constants.h"
#include "Indexing.h"
#include "Coordinates.h"

#include "Scatterer.h"

// Super-class for incoming waves
class IncWave
{
public:
static Indexing global_index;

enum Type{PW, Pt};      // type of wave

std::vector<complex> coeff;
complex k_out;

// Constructor/Destructor
IncWave(){};
~IncWave(){};
  
  // Translation of expansion to a different origin
 virtual void Translate( Pvec& c, std::vector<complex>& T_c) = 0;
 //virtual Pvec getLocation() = 0;
 virtual complex Evaluate(Pvec& p) = 0;
};



// Spherical wave
class SphericalWave : public IncWave{

Indexing* SW_index;
complex strength;
Pvec loc;

public:
SphericalWave(complex k_out_, complex strength_, const Pvec& loc_);

~SphericalWave(){ delete SW_index; }
      
 void Translate( Pvec& c, std::vector<complex>& T_c);
 complex Evaluate(Pvec& p);
 Pvec getLocation() { return loc; }
 complex getAmplitude() { return strength; }
  
};



// Plane wave
class PlaneWave : public IncWave{
  
 public:
  complex amp;     // amplitude
  Pvec wv;         // wave vector
  
  PlaneWave(complex amp_, Pvec wv_);
  ~PlaneWave(){}
  
  void Translate( Pvec& c, std::vector<complex>& T_c);
  complex Evaluate(Pvec& p);
  Pvec getWaveVector() { return wv; }
  complex getAmplitude() { return amp; }
};



#endif
