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


//Indexing PSC_Env::global_index(LMAX);

void PSC_Env::Construct( complex k_out_, std::vector<Scatterer> ScL_ )
{
  // Wave numbers
  k_out = k_out_;  

  // TODO: streamline. Really needed?
  // Change format for the list of scatterers
  ScL.resize(ScL_.size());
  for( int i = 0; i < (int) ScL_.size(); i++ )
    ScL[i].copy(ScL_[i]);
  
  // Build transfer operator
  transfer = new Transfer(nLevels, eps, k_out, ScL);
}

// Safe releasing of memory
void PSC_Env::Destruct(){
  FMM::S_Rotation::S_Delete();
  FMM::Z_transfer::Z_Delete();

  delete transfer; transfer = NULL;
}








#endif
