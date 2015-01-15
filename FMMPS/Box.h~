#ifndef BOX_FMM_H
#define BOX_FMM_H

#include "General.h"
#include "Vec3.h"
//#include "NFunction.h"

// *** Only contains .h file ***

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
