/*
 *  general.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  All libraries and functions required by most files
 */

#ifndef INDEXING_H
#define INDEXING_H


//---------Indexing----------//

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
