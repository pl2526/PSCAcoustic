/*
 *  Solver.h
 *  PSCAcoustic
 *
 *  Iterative solver using Krylov subspace method
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

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <complex>
#include <iostream>

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "PSCEnv.h"


#include "./Aux/Eigen/Dense"

using namespace Eigen;


class Solver
{

Indexing Index;

public:

// Constructor/Destructor
Solver();
~Solver(){}

// Initialization of the right-hand side
void RHSInit(IncWave* IW, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& RHS);

 void RHSInit(std::vector<IncWave*> IW, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& RHS);

// Fast application of the operator TA
 void ScatteringOperator(PSC_Env& PSC_env, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& u,
			 std::vector< std::vector<complex > >& v);

// Solution of the linear system
 void Solve(PSC_Env& PSC_env, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& RHS, 
	    std::vector< std::vector<complex> >& u,  double& res, double& rel_res, double& cond, int& Niter, int idx = 0  );

};

#endif
