#ifndef PSCEnv_CPP
#define PSCEnv_CPP

/*
 *  PSCEnv.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Set up the environment for the Particle Scattering (ADD DETAILS)
 */

#include "PSCEnv.h"


// TODO: Merge Transfer and PSC env
// TODO: nLevels should be passed

void PSC_Env::Construct( int nlevels_, double r, complex k_in, complex k_out_, std::vector<Scatterer> ScL_ )
{
  nlevels = nlevels_;
  k_out = k_out_;
  cout << "nlevels : " << nlevels << endl;


  // TODO: streamline. Really needed?
  // Change format for the list of scatterers
  ScL.resize(ScL_.size());
  for( int i = 0; i < (int) ScL_.size(); i++ )
    ScL[i].copy(ScL_[i]);
  
  // Build transfer operator
  transfer = new Transfer(nlevels, eps, r, k_in, k_out, ScL);
}

// Safe releasing of memory
void PSC_Env::Destruct(){
  FMM::S_Rotation::S_Delete();
  FMM::Z_transfer::Z_Delete();

  delete transfer; transfer = NULL;
}








#endif
