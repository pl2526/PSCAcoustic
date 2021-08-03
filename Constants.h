/*
 *  Constants.h
 *  PSCAcoustic
 *
 *  All sorts of static constants.
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


#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>

typedef std::complex<double> complex;

//Constants
const static int DIM = 3;
const double  PI(M_PI);
const complex CI(0,1);

//-----------Parameters------------//

//Convention for spherical coordinates
// r >= 0
// 0 <= theta <= PI (inclination)
// 0 <= phi <= 2*PI (azimuth)

// Output filename
const std::string LOC_PREFIX("_Loc_");
const std::string INFO_PREFIX("_Info_");
const std::string SOL_PREFIX("_Sol_");

// Flags
const bool CT_PRECOMPUTE = true;   // Pre-compute close-term transfer functions
const bool FAST_TRANSFER = true;   // Perform transfer using fast algorithm or not

// Display switch
// 0: no display
// 1: timing
// 2: Initialization progress, timing
const int verbose = 2;

// KSP solver (PetSc)
const double RTOL = 1e-5;
const double ATOL = 1e-14;
const double DTOL = 1e5;
const int MAXITS = 650;
const int K_LS = 20;

// Background properties
const double OMEGA = 1.;   // Angular frequency (Hz)

const double RHO_OUT = 1.; // Density 
const double C_OUT = 1.;    // Longitudinal velocity
const complex K_OUT = 2*PI*OMEGA/C_OUT;    // Longitudinal wave number


// Scatterers properties
const int NScat = 2000;     // Total number of scatterers (set in PSCEnv.h)

const double RADIUS = 0.0053;
const double D = 1.;               // Diameter of computational domain
const double D_MIN = RADIUS * 0.5;

// *** CURRENTLY NO LAME CTS for scatterers; assumed rigid

const double RHO = 0.2; // Density 
const double C = 0.1;    // Longitudinal velocity
const complex K = 2*PI*OMEGA/C;    // Longitudinal wave number



// Parameters for FMM
const int nLevels = 3;      // Number of levels in FMM (Should be >= 3)
const double EPS = 1e-6;          // Desired accuracy for FMM


const double CUTOFF_SIZE = 2. * C_OUT/OMEGA;   // Cutoff boxsize (in number of wavelengths) for interface between low- and high-frequency regimes



// TODO: SHOULD BE OBTAINED FROM NUMERICAL ANALYSIS
const int LMAX = 1;        // Maximum degree for representation of Cartesian coordinates in Spherical Wave Functions


enum SWF_Type{Hankel, Bessel};

#endif
