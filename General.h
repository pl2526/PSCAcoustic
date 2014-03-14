/*
 *  general.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  All libraries and functions required by most files
 */

#ifndef GENERAL_H
#define GENERAL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <assert.h>
#include <time.h>

#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <complex>  


// Armadillo linear algebra library
//#include "../Armadillo/usr/include/armadillo"
#include "../armadillo-4.000.2/include/armadillo"

typedef std::complex<double> complex;
typedef std::ostream ostream;
typedef std::fstream fstream;
typedef std::ofstream ofstream;
typedef std::ios ios;
#define cout std::cout
#define cerr std::cerr
#define endl std::endl
#define map std::map

#include "Constants.h"
#include "SpecialFunctions.h"



#endif
