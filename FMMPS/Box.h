/*
 *  Coordinates.cpp
 *  PSCAcoustic
 *
 *  Box of the FMM tree.
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




#ifndef BOX_FMM_H
#define BOX_FMM_H

#include "General.h"
#include "Vec3.h"

namespace FMM {
  
  class Box 
  {  
    std::vector<complex> M;             // Multipole
    std::vector<complex> L;             // Local expansion
    
    
  public:
    typedef std::vector<int> PointList;
    typedef PointList::iterator pointIter;
    
    int n;                   // The box number
    Vec3 center;             // The box center
    
    Box* parent;             // A pointer to this box's parent
    PointList pointIndex;    // If leaf, contains indices of points it contains
    
    Box() {}
  Box( int n_, const Vec3& center_ ) : n(n_), center(center_) {}
    ~Box() {}
    
    inline void addPointIndex( int index ) {
      pointIndex.push_back( index );
    }
    
    //This should set the size of the expansion only
    inline void makeMultipole( int size ) {
      M.resize( size );
    }
    
    inline void makeLocal( int size ) {
      L.resize( size );
    }
    
    inline std::vector<complex>& getMultipole() {
      return M;
    }
    
    inline std::vector<complex>& getLocal() {
      return L;
    }
    
    
  private:
    // Disable Copy and Assignment
    Box( const Box& S ) {}
    void operator=( const Box& S ) {}
  };
  
}

#endif
