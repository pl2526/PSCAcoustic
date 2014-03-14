#ifndef TRANSFER_CPP
#define TRANSFER_CPP

/*
 *  Transfer.cpp
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Implementation of the transfer pass; to transfer a given
 *  field expansion (in spherical harmonics) from one center 
 *  to a different center of expansion
 *
 */

#include "Transfer.h"


//Constructor
Transfer::Transfer( int nLevels, double eps, complex  k_out_, std::vector<Scatterer>& ScL) : k_out(k_out_)
{
  
  // set wave number for Helmholtz kernel
  KERNEL.kappa = k_out;
  
  p.resize(NScat);
  for( int i = 0; i < NScat; i++ ){
    FMM::Vec3 FMMVec( ScL[i].getLocation().x, ScL[i].getLocation().y, ScL[i].getLocation().z );
    p[i] =  FMMVec;
  }
  
  //Build FMM environment (only for multiple scatterers)
  if( FAST_TRANSFER && NScat > 1 ){
    MLFMMEnv = new FMM::MLFMM_Env(KERNEL, k_out, p, ScL, nLevels, eps);
  } else {
    //TODO : Need to take inc wave into consideration...
    FMM::S_Rotation::S_Init( 50 );
  }
  
}


//Destructor
Transfer::~Transfer(){
  if( FAST_TRANSFER && NScat > 1 ){
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
  if( NScat > 1 )
    MLFMMEnv->execute( PWE, T_PWE, k_out, FMM::FORWARD);
}


#endif


















