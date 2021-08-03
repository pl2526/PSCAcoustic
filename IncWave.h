/*
 *  IncWave.h
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
