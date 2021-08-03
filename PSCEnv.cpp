#ifndef PSCEnv_CPP
#define PSCEnv_CPP

/*
 *  PSCEnv.h
 *  PSC
 *
 *  Set up the environment for the Particle Scattering
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
