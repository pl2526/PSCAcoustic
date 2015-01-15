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
