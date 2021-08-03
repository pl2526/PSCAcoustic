/*
 *  Transfer.h
 *  PSCAcoustic
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


















