/*
 *  PSCEnv.h
 *  PSCAcoustic
 *
 *  Set up the environment for fast matrix-vector product in the T-matrix context.
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


#ifndef PSCEnv_H
#define PSCEnv_H

#include "General.h"
#include "Scatterer.h"
#include "Transfer.h"
#include "./FMMPS/TransferUtils.h"


class PSC_Env
{

 public:

  //static Indexing global_index;

  std::vector<Scatterer> ScL;       // Scatterer list
  Transfer *transfer;               // Longitudinal waves

  int nlevels;
  double eps;
  complex k_out;

  // Constructor
 PSC_Env(double eps_) : eps(eps_) {}

  // Empty destructor 
  ~PSC_Env(){}


  //Accessors
  inline std::vector<Scatterer> getScatList() { return ScL; }
  inline Transfer* getTransfer() { return transfer; }

  // Initializer
  std::vector<Scatterer> ScatInit( int N  = 1, complex k = 2*PI );   //Initialize location and properties of scatterers

 public:
  void Construct(int nlevels, double r, complex k_in, complex  k_out_, std::vector<Scatterer> ScL_ );    // Initialization of environment
  void Destruct();                                                                        // Clear environment
  void Finalize( std::vector< std::vector<complex > >& sol_out,
		 std::vector< std::vector<complex > >& sol_in, Pvec& source_loc, 
		 int nIter, double res_norm, double rel_res_rhs, int nProc);              // Post-processing of solution

};



#endif
