/*
 *  General.h
 *  PSCAcoustic
 *
 *  Libraries and functions required by most files
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



#ifndef GENERAL_H
#define GENERAL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <assert.h>
#include <time.h>
#include <mpi.h>

#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <complex>  


// Armadillo linear algebra library
//#include "../Armadillo/usr/include/armadillo"
#include "./armadillo-4.000.2/include/armadillo"

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
