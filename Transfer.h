/*
 *  Transfer.h
 *  PSCAcoustic
 *
 *  Implementation of the transfer pass; to transfer a given
 *  field expansion (in spherical harmonics) from one center 
 *  to a different center of expansion
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
  Transfer( int nLevels, double eps, double r, complex k_in, complex  k_out, std::vector<Scatterer>& ScL);

  //Destructor
  ~Transfer();

  //Execution of transfers
  void execute( const std::vector< std::vector<complex > >& PWE, std::vector< std::vector<complex > >& T_PWE );
};


#endif


















