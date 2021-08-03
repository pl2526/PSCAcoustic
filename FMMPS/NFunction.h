/*
 *  NFunction.cpp
 *  PSCAcoustic
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

#ifndef NFUNCTION_FMM_H
#define NFUNCTION_FMM_H

#include "General.h"
#include "Quadrature.h"


namespace FMM {
  
  // Super class for scratch space
  //TODO :could it be a struct?
  class NFunction
  {
  public:
    
    std::vector<complex> C;
    
    // This = 0
    inline void zero()
    {
      C.assign( C.size(), 0 );
    }

    inline void setSize( int size )
    {
      C.resize( size );
    }

    inline int size()
    {
      return C.size();
    }

    inline complex operator [](int i) { return C[i]; }

  };




  //TODO : this is HF scratch space. Kept old name for convenience. Must change later
  class HF_NFunction : public NFunction 
    {
    public:
      
      Quadrature* quad;     // Pointer to quadrature of field values
      
      // Construct a numerical function over a given quadrature
    HF_NFunction( Quadrature* q_ = NULL ) : NFunction() { setQuad(q_); }
      // Destructor
      ~HF_NFunction()
	{
	  // Don't delete quad, it doesn't belong to us
	}
      
    //Set quadrature and size of std::vectors holding expansions
      inline void setQuad( Quadrature* q )
      {
	quad = q;
	C = std::vector<complex>( quad == NULL ? 0 : quad->size() );
      }
      
      inline Quadrature& getQuad() 
      { 
	return *quad;
      }
    };
  
}

#endif
