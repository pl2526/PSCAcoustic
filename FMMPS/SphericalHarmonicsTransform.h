/*
 *  SphericalHarmonicsTransform.h
 *  PSCAcoustic
 *
 *  Fast spherical harmonics transform.
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

#ifndef SPHERICAL_HARMONICS_TRANSFORM_H
#define SPHERICAL_HARMONICS_TRANSFORM_H

#include "../Indexing.h"

namespace FMM {
  
  class SH_Transform
  {
    Indexing* index;
    
    std::vector< std::vector<complex> > b;
    std::vector< std::vector<complex> > c;
    std::vector<complex> vec1;
    complex* array;
    std::vector< std::vector< std::vector<complex> > > Legendre;
    std::vector<double> theta;
    std::vector< std::vector<double> > weight;
    std::vector< std::vector<double> > phi;
    std::vector<int> M;
    int L, Q, N;

    
  public:

    // Constructor/Destructor
    SH_Transform(Quadrature* quad, Indexing* index);
    ~SH_Transform();
    
    // Transforms    
    std::vector<complex> Forward( std::vector<complex>& vec2 ); // Forward transform
    std::vector<complex> Backward( std::vector<complex>& vec2 ); // Backward transform
    
  };
  

}

#endif
