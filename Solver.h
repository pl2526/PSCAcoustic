#ifndef SOLVER_H
#define SOLVER_H

/*
 *  Solver.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Iterative linear solver implementation 
 *  through PETSc library.
 */

#include <vector>
#include <complex>
#include <iostream>
//#include <mpi.h>

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "PSCEnv.h"


#include "./Eigen/Dense"

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
	    std::vector< std::vector<complex> >& u,  double& res, double& rel_res, double& cond );

};

#endif