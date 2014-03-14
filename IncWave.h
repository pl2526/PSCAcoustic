#ifndef INCWAVE_H
#define INCWAVE_H

/*
 *  IncWave.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Incoming wave multipole expansion and satellite functions
 */

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
