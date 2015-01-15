/*
 *  NFunction.cpp
 *  PSCAcoustic
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
