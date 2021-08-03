/*
 *  PSCEnv.h
 *  PSCAcoustic
 *
 *  Set up the environment for fast matrix-vector product in the T-matrix context.
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
