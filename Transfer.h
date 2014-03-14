/*
 *  Transfer.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Implementation of the transfer operator
 *  for two different speed of propagation.
 *
 */

#ifndef TRANSFER_H
#define TRANSFER_H

#include "./FMMPS/General.h"
#include "./FMMPS/MLFMM_Env.h"

// Transfer is oblivious to the fact that the speed of sound is different for different kinds of waves
  
class Transfer
{

  //Objects needed by FMM
  FMM::HelmKernel KERNEL;    // Longitudinal waves
  FMM::MLFMM_Env *MLFMMEnv;

  std::vector<FMM::Vec3> p;
  complex  k_out;

 public:
  
  //Constructor
  Transfer( int nLevels, double eps, complex  k_out, std::vector<Scatterer>& ScL);

  //Destructor
  ~Transfer();

  //Execution of transfers
  void execute( const std::vector< std::vector<complex > >& PWE, std::vector< std::vector<complex > >& T_PWE );
};


#endif


















