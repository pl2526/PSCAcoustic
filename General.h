/*
 *  General.h
 *  PSCAcoustic
 *
 *  Libraries and functions required by most files
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
