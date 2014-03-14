#ifndef PSCEnv_H
#define PSCEnv_H

/*
 *  PSCEnv.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Set up the environment for the Particle Scattering (ADD DETAILS)
 */

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

  int nLevels;
  double eps;
  complex k_out;

  // Constructor
 PSC_Env( int nLevels_, double eps_) : nLevels(nLevels_), eps(eps_) {}

  // Empty destructor 
  ~PSC_Env(){}


  //Accessors
  inline std::vector<Scatterer> getScatList() { return ScL; }
  inline Transfer* getTransfer() { return transfer; }

  // Initializer
  std::vector<Scatterer> ScatInit( int N  = 1, complex k = 2*PI );   //Initialize location and properties of scatterers

 public:
  void Construct( complex  k_out_, std::vector<Scatterer> ScL_ );    // Initialization of environment
  void Destruct();                                                                        // Clear environment
  void Finalize( std::vector< std::vector<complex > >& sol_out,
		 std::vector< std::vector<complex > >& sol_in, Pvec& source_loc, 
		 int nIter, double res_norm, double rel_res_rhs, int nProc);              // Post-processing of solution

};



#endif
