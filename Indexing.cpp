/*
 *  general.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  All libraries and functions required by most files
 */

#ifndef INDEXING_CPP
#define INDEXING_CPP


#include <vector>
#include <cmath>
#include <complex>  

typedef std::complex<double> complex;

#include "Indexing.h"


// TODO: Improve efficiency

// Constructors
Indexing::Indexing( int lmax ) { Compute(lmax); }

void Indexing::Compute( int lmax ) { 
  int Max =  (lmax+2)*(lmax) + 1; //Total number of harmonics of order <=lstd::max
  //std::vector< std::vector<int> > v(Max, std::vector<int>(2) );
  Index.resize(Max);
  
  int i;
  for ( int l=0 ; l<=lmax ; l++){
    for ( int m=-l ; m<=l ; m++){
      i = std::max((l+1)*(l-1)+1 , 0) + (m+l);

      Index[i].resize(2);
      Index[i][0] = l;
      Index[i][1] = m;
      //v[i][0] = l;
      //v[i][1] = m;
    }
  }
  
  //Index = v;
}


#endif
