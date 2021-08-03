/*
 *  general.h
 *  PSC
 *
 *  All libraries and functions required by most files
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
