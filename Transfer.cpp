/*
 *  Transfer.cpp
 *  PSC
 *
 *  Implementation of the transfer pass; to transfer a given
 *  field expansion (in spherical harmonics) from one center
 *  to a different center of expansion
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


#ifndef TRANSFER_CPP
#define TRANSFER_CPP

#include "Transfer.h"


//Constructor
Transfer::Transfer( int nlevels, double eps, double r, complex k_in, complex  k_out_, std::vector<Scatterer>& ScL) : k_out(k_out_)
{
  
  // set wave number for Helmholtz kernel
  KERNEL.kappa = k_out;
  
  p.resize(ScL.size());
  for( int i = 0; i < (int) ScL.size(); i++ ){
    FMM::Vec3 FMMVec( ScL[i].getLocation().x, ScL[i].getLocation().y, ScL[i].getLocation().z );
    p[i] =  FMMVec;
  }
  
 cout << "foo1" << endl;

  //Build FMM environment (only for multiple scatterers)
  if( FAST_TRANSFER && p.size() > 1 ){
    MLFMMEnv = new FMM::MLFMM_Env(KERNEL, r, k_in, k_out, p, ScL, nlevels, eps);
  } else {
    //TODO : Need to take inc wave into consideration...
    FMM::S_Rotation::S_Init( 50 );
  }
  
}


//Destructor
Transfer::~Transfer(){
  if( FAST_TRANSFER && p.size() > 1 ){
    delete MLFMMEnv;
    MLFMMEnv = NULL;
  } else {
    FMM::S_Rotation::S_Delete();
  }
}


//Execute transfer of multipole expansions
// void execute( const complex  outgoing[], complex  transfered[], FMM::Type type)
void Transfer::execute( const std::vector< std::vector<complex > >& PWE, std::vector< std::vector<complex > >& T_PWE)
{
  if( p.size() > 1 )
    MLFMMEnv->execute( PWE, T_PWE, k_out, FMM::FORWARD);
}


#endif


















