/*
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

#ifndef TRANSFER_VECTOR_FMM_H
#define TRANSFER_VECTOR_FMM_H

#include "General.h"
#include "Box.h"
#include "Transfer_Function.h"

// ***Only contains .h file

namespace FMM {
  
  struct Close_Vector 
  {
    Box* b1;
    Box* b2;
    
    Close_Vector( Box* b1_, Box* b2_ ) : b1(b1_), b2(b2_) {}
    
    friend ostream& operator<<(ostream& os, const Close_Vector& t)
    {
      return (os << t.b1->n << "-" << t.b2->n);
    }
  };
  
  
  struct Trans_Vector 
  {
    Box* b1;                  // Pointer to box 1 (from box)
    Box* b2;                  // Pointer to box 2 (to box)
    
    Transfer_Function* T;     // Pointer to Transfer Function to be applied
    
    int x, y, z;              // Representation of r_0 in number of boxes
    
  Trans_Vector() : b1(NULL), b2(NULL), T(NULL), x(0), y(0), z(0) {}
    
  Trans_Vector( Box* b1_, Box* b2_, Vec3 r0_ )
  : b1(b1_), b2(b2_), T(NULL), x(r0_.x), y(r0_.y), z(r0_.z) {}
    
    // Comparator for sorting
    inline bool operator<(const Trans_Vector& t) const {
      return (x < t.x 
	      || (x == t.x && (y < t.y
			       || (y == t.y && z < t.z))));
    }
    
    friend ostream& operator<<(ostream& os, const Trans_Vector& t)
    {
      return (os << t.b2->n << "->" << t.b1->n << ":"
	      << "\t(" << t.x << ", " << t.y << ", " << t.z << ")");
    }
  };
  
}

#endif

