/*
 *  Finalize.h
 *  PSCAcoustic
 *
 *  Post-processing for writing solution to file in 
 *  standard form.
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


#ifndef FINALIZE_H
#define FINALIZE_H

#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "Scatterer.h"
#include "IncWave.h"

// Write to file information pertaining to the problem solved
void Write_Info(int Proc_Idx, std::string& type,double residual, double relative_residual, double cond, int Niter, std::string filename, complex k_out = K_OUT, complex k_in = K, int N = NScat);

// Write to file information about each source
void Write_Source( int Proc_Idx, IncWave* IW, IncWave::Type type, std::string filename);
void Write_Source( int Proc_Idx, std::vector<IncWave*>& IW, IncWave::Type type, std::string filename);

// Write to file location of each scatterer
void Write_Location( int Proc_Idx, std::vector<Scatterer>& ScL, std::string filename);

// Write to file solution in VSWF basis
void Write_Solution( int Proc_Idx, Indexing& index, std::vector< std::vector<complex> > u, std::string filename);


#endif
