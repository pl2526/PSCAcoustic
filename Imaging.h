/*
 *  Imaging.h
 *  PSCAcoustic
 *
 *  Produces image of wave field.
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


#ifndef IMAGING_H
#define IMAGING_H

#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"

// u_in: RSWF expansion for field inside scatterers

void Imaging( std::vector<Scatterer>& ScL, IncWave* IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL, complex k_out = K_OUT, 
	      complex k_in = K);

void Imaging( std::vector<Scatterer>& ScL, std::vector<IncWave*>& IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL, complex k_out = K_OUT, 
	      complex k_in = K);


#endif
