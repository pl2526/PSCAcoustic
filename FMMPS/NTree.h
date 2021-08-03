/*
 *  NTree.h
 *  PSCAcoustic
 *
 *  FMM Tree.
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


#ifndef NTREE_FMM_H
#define NTREE_FMM_H

#include "General.h"
#include "Vec3.h"
#include "Box.h"

namespace FMM {
  
  // A class templated to create
  // DIM = 1: BinaryTrees    Vec3's are (x,0,0)
  // DIM = 2: QuadTrees      Vec3's are (x,y,0)
  // DIM = 3: OctTrees       Vec3's are (x,y,z)
  
  
  //template <int DIM>
  class NTree
  {
    double rootSize;             // The size of the root box
    Vec3 pmin;                   // Handle to the bottom left of the root box
    
    int nLevels;                 // The number of levels: 0 = root, nLevel = leaf
    int nBoxes;
    
    // Should sparsify this if very few boxes are actually required
    std::vector<Box*> boxes;          // The boxes in this tree
    
  public:
    
    const static int BRANCH = 1 << DIM;   // The branching factor = 2^DIM  

    // Constructor of an NTree with levels [0,L]
    NTree( std::vector<Vec3>& p, int L );
    
    
    // Destructor
    ~NTree();
    
    
    // ---Accessors--- //
    
    inline int numLevels(){ return nLevels; }   
    inline Box* getBox( int n, int L ){ return boxes[ getBoxIndex(n,L) ]; }
    
    // Note that L = 0 is the root
    inline double getBoxSize( int L ){ return rootSize / (1 << L); } 
    
    // The kth-order parent of box n
    inline int parent( int n, int k = 1 ){ return n >> (DIM*k); }
    
    // The mth kth-order child of box n
    inline int child( int n, int m = 0, int k = 1 ){ assert( m < (1 << (DIM*k)) ); return m + (n << (DIM*k)); }
    
    int getBoxIndex( int n, int L );
    Vec3 boxCenter( int n, int L );
    
    
    
  private:
    int interleave( int x, int y, int z );
    Vec3 deinterleave( int n );  
    std::vector<int> siblings( int n );

    // Disable Copy and Assignment
    NTree( const NTree& S ) {}
    void operator=( const NTree& S ) {}



  public:

    // ---Output--- //
    
    friend ostream& operator<<(ostream& os, NTree& t)
    {
      int nLevels = t.numLevels();
      for( int n = 0; n < (1 << (DIM*nLevels)); ++n ) {
	Box* b = t.getBox(n,nLevels);
	if( b == NULL ) continue;
	
	for( int l = nLevels; l >= 0; --l ) {
	  os << t.parent(b->n,l) << "\t";
	}
	os << b->center << endl;      
      }
      return os;
    }
    
  };
  
}

#endif
