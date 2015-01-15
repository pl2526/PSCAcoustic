/*
 *  Distributions.h
 *  PSCAcoustic
 *
 *  Various routines to generate arrangements of scatterers.
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


#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "Constants.h"
#include "Indexing.h"
#include "Scatterer.h"

// TODO: Comment

// Prototype for arrangement generator
std::vector<Scatterer> ScatInit( int N, complex k, complex k_out );

// Prototype for various distributions
std::vector<Pvec> ConcentricSpheres(double R, double r, double d);
std::vector<Pvec> CylSlab(Pvec a, double r, double R, double h);
std::vector<Pvec> RandSpherical(double R, double r, double d, Pvec c, double N);
std::vector<Pvec> TruncatedSphere(double R, double t, double r, double d, Pvec c, double N, Pvec src_loc, bool rand_flag = false);
std::vector<Pvec> RandSphericalXY(double R, double r, double d, Pvec c, double N);
std::vector<Pvec> RandRectangularXY(double x_dim, double y_dim, double r, double d, Pvec c, double N);
std::vector<Pvec> RandCubic(double Dx, double Dy, double Dz, double r, double d, Pvec c, int N);
std::vector<Pvec> SphericalPeriodic(double R, Pvec center, double r, double p);
std::vector<Pvec> ErgodicCavity(double R, double t, double thick, Pvec center, double r, double d, double N, Pvec src_loc, bool rand_flag);
std::vector<Pvec> ErgodCavityPeriodic(double R, double t, Pvec center, double r, double p);
std::vector<Pvec> ErgodCavity2(double R, Pvec center, double r, double d, double N);




#endif
