/*
 *  Coordinates.cpp
 *  PSCAcoustic
 *
 *  Box of the FMM tree.
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
