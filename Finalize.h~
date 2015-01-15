#ifndef FINALIZE_H
#define FINALIZE_H

/*
 *  Finalize.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 10/3/13.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  All post-processing for writing solution to file in 
 *  standard form.
 */


#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "Scatterer.h"
#include "IncWave.h"

// Write to file information pertaining to the problem solved
void Write_Info(int Proc_Idx, std::string& type,double residual, double relative_residual, double cond, std::string filename);

// Write to file information about each source
void Write_Source( int Proc_Idx, IncWave* IW, IncWave::Type type, std::string filename);
void Write_Source( int Proc_Idx, std::vector<IncWave*>& IW, IncWave::Type type, std::string filename);

// Write to file location of each scatterer
void Write_Location( int Proc_Idx, std::vector<Scatterer>& ScL, std::string filename);

// Write to file solution in VSWF basis
void Write_Solution( int Proc_Idx, Indexing& index, std::vector< std::vector<complex> > u, std::string filename);


#endif
