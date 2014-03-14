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
const double RTOL = 1e-11;
const double ATOL = 1e-13;
const double DTOL = 1e5;
const int MAXITS = 200;
const int K_LS = 5;

// Background properties
const double OMEGA = 1.;   // Angular frequency (Hz)

const double RHO_OUT = 1.; // Density 
const double C_OUT = 1.;    // Longitudinal velocity
const complex K_OUT = 2*PI*OMEGA/C_OUT;    // Longitudinal wave number


// Scatterers properties
const int NScat = 1000;     // Total number of scatterers (set in PSCEnv.h)

const double RADIUS = 0.00001;//0.4;
const double D = 1.;               // Diameter of computational domain
const double D_MIN = RADIUS * 0.5;

// *** CURRENTLY NO LAME CTS for scatterers; assumed rigid

const double RHO = 5;//1.05;//5.; // Density 
const double C = 0.11;//0.85;//0.11;    // Longitudinal velocity
const complex K = 2*PI*OMEGA/C;    // Longitudinal wave number



// Parameters for FMM
const int nLevels = 3;      // Number of levels in FMM (Should be >= 3)
const double EPS = 1e-12;          // Desired accuracy for FMM


const double CUTOFF_SIZE = 2. * C_OUT/OMEGA;   // Cutoff boxsize (in number of wavelengths) for interface between low- and high-frequency regimes



// TODO: SHOULD BE OBTAINED FROM NUMERICAL ANALYSIS
const int LMAX = 0;//6;//14;//10;        // Maximum degree for representation of Cartesian coordinates in Spherical Wave Functions


enum SWF_Type{Hankel, Bessel};

#endif
