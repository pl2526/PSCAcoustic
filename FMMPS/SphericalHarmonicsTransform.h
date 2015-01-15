/*
 *  SphericalHarmonicsTransform.h
 *  PSCAcoustic
 *
 *  Fast spherical harmonics transform.
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
