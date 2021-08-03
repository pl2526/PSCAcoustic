

#ifndef SR_CPP
#define SR_CPP

/*
 *  MAIN.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Overstructure.
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

#include <vector>
#include <complex>
#include <iostream>

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "PSCEnv.h"
#include "Solver.h"
#include "Imaging.h"
#include "Finalize.h"

Indexing temp_index(LMAX);



int main(int argc,char **args)
{
  int Proc_Idx = 0;


  // Properties of source
  std::string src_type("PtW");
  //double amplitude = 1;
  Pvec src_location(25., 30., 40.,Pvec::Cartesian);
  complex src_amplitude = 1.;


  // Properties of images
  double L = 1.5;    // Size of window [-L,L]x[-L,L]
  int M = 65;         // Number of samples per dimension (2*M+1)

  // Display relevant information concerning the problem
  cout << endl << endl;
  cout << "LMAX : " << LMAX << endl;
  cout << "RTOL : " << RTOL << endl;
  cout << "MAXITS : " << MAXITS << endl;
  cout << "OMEGA : " << OMEGA << endl;
  cout << "C_OUT : " << C_OUT << endl;
  cout << "C : " << C << endl;
  cout << "K_OUT : " << K_OUT << endl;
  cout << "K : " << K << endl;
  cout << "NScat : " << NScat << endl;
  cout << "RADIUS : " << RADIUS << endl;
  cout << "nLevels : " << nLevels << endl;
  cout << "EPS : " << EPS << endl;

  cout << endl;

  // PSC environment for fast algorithm
  PSC_Env PSC_env(nLevels, EPS); 

  // Linear solver
  Solver solver; 

  // TODO: Should not have different instances of Index
  Indexing Index(LMAX);

  // Construct scatterers
  cout << "***Building scatterers..." << endl;
  std::vector<Scatterer> ScL = ScatInit(NScat, K, K_OUT );
  cout << "***Building scatterers: done" << endl << endl;

  // Construct PSC environment
  cout << "***Constructing environment..." << endl;
  PSC_env.Construct(K_OUT, ScL);
  cout << "***Constructing environment: done" << endl << endl;
 
  



  // *** Solve forward problem

  // Initialization of source
  // R.h.s.
  SphericalWave* IW = new SphericalWave(K_OUT, src_amplitude, src_location);
  std::vector< std::vector<complex > > RHS(NScat, std::vector<complex >(Index.size())); 
  solver.RHSInit(IW, ScL, RHS);  

  // Solve linear system
  cout << "   ***Solving linear system..." << endl;
  std::vector< std::vector<complex> > u(NScat, std::vector<complex>(Index.size()));
  double res = 1e10, rel_res = 1e10, cond;
  solver.Solve(PSC_env, ScL, RHS, u, res, rel_res, cond);
cout << "   ***Solving linear system: done" << endl;


// Compute solution inside scatterers
std::vector<std::vector<complex> > *u_in = new std::vector<std::vector<complex> >(NScat, std::vector<complex>((LMAX+1)*(LMAX+1)));

  // Produce image
  std::string im_filename("SuperRes_Im");
  Imaging(ScL, IW, im_filename, u, Proc_Idx, M, L, u_in);
  
  // Write information about problem to file
  cout << "   ***Writing to file..." << endl;
  std::string filename("Super_Res_A");
  Write_Info(Proc_Idx, src_type, res, rel_res, cond, filename);
  Write_Source(Proc_Idx, IW, IncWave::Pt, filename);
  Write_Location(Proc_Idx, ScL, filename);
  Write_Solution( Proc_Idx, Index, u, filename);
  cout << "   ***Writing to file: done" << endl << endl;
  


  





  return 0;
}






#endif
