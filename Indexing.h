/*
 *  Indexing.h
 *  PSCAcoustic
 *
 *  Linear indexing for spherical wave functions and spherical harmonics.
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

#ifndef INDEXING_H
#define INDEXING_H


class Indexing
{
  std::vector< std::vector<int> > Index;

 public:
  // Constructors
  Indexing(){}
  Indexing( int lmax );

  //Empty destructor
  ~Indexing(){}

  // Initialization
  void Compute( int lmax );
 
  // TODO: clean up from here

  // Accessors
  inline std::vector< std::vector<int> > getIndex(){ return Index; }
  inline int LocateIndex(int l, int m){ return l*l + (l+m); }
  static inline int Idx(int l, int m){ return l*l + (l+m); }
  inline int size() const { return Index.size(); }
  
  // Operators (overloaded)
  inline std::vector<int> operator ()(int l) { return Index[l]; }
  inline int operator ()(int l, int m) { return Index[l][m]; } 
  inline void operator =(const Indexing& idx) { 
    Index.resize(idx.Index.size());
    for( int i = 0; i < (int) idx.Index.size(); i++ ){
      Index[i][0] = idx.Index[i][0];
      Index[i][1] = idx.Index[i][1];
    }
  }

};




#endif
