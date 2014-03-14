#ifndef IMAGING_H
#define IMAGING_H

/*
 *  Imaging.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 11/11/13.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Produce image of wave field
 */



#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"

// u_in: RSWF expansion for field inside scatterers

void Imaging( std::vector<Scatterer>& ScL, IncWave* IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL);

void Imaging( std::vector<Scatterer>& ScL, std::vector<IncWave*>& IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL);


#endif
