/*
 *  Indexing.h
 *  PSCAcoustic
 *
 *  Linear indexing for spherical wave functions and spherical harmonics.
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
