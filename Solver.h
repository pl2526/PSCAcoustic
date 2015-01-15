/*
 *  Solver.h
 *  PSCAcoustic
 *
 *  Iterative solver using Krylov subspace method
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
