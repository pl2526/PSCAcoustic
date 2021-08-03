/*
 *  Distributions.h
 *  PSCAcoustic
 *
 *  Various routines to generate arrangements of scatterers.
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
