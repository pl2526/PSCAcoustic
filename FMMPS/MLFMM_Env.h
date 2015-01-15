/*
 *  MLFMM_Env.h
 *  PSCAcoustic
 *
 *  Multi-level FMM environment
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


#ifndef MLFMM_Env_FMM_H
#define MLFMM_Env_FMM_H

#include "General.h"
#include "NTree.h"
#include "Level.h"
#include "SphericalHarmonicsTransform.h"
#include "../Scatterer.h"
#include "../Indexing.h"


namespace FMM {
  
  // Determines interface between low- and high-frequency regimes
  inline bool Interface( double boxSize, double cutoff ){ return ( boxSize <= cutoff && 2*boxSize > cutoff ) ? true : false; }
  


  class MLFMM_Env
  {
  public:
    static Indexing global_index;
    HelmKernel K;
    complex k_out;
    complex k_in;
    std::vector<Scatterer> ScL;
    
    std::vector<Leaf_Translation_Function*> leaf_Translation;
    std::vector<Vec3> sourcefield;
    std::vector< std::vector<complex> > SfH;
    SH_Transform* leaf_SH_T;
    
    NTree tree;
    std::vector<Level> levels;
    
    // The interface level is defined as the first level in the HF regime with box size above the cutoff
    // or the leaf level if it is already in the HF regime
    int Interface_Level;
    
    
    // Constructor
    // TODO : Should pass ScL by reference but need to be constant...
    MLFMM_Env( HelmKernel K, double r, complex k_in, complex k_out, std::vector<Vec3> p, std::vector<Scatterer> ScL, int levels_, double eps);
    
    // Destructor
    ~MLFMM_Env();
    
    // Computes   omega = A*psi
    void execute( const std::vector< std::vector<complex> >& PWE, std::vector< std::vector<complex> >& T_PWE, complex kappa, Type type);
    
    
    
    
    //-------- Private --------//
    
    // Should be private...
    inline Level& getLevel( int L ) { return levels[L]; }
    inline Box* getBox( int n, int L ) { return tree.getBox(n,L); }
    
  private:
    
    // Main algorithm to determine transfer pairs and close pairs
    // Two boxes are transfers if their parents are neighbors and they aren't
    // Two boxes are close if they are leaf level and are neighbors  
    void defineTransferAndClose( Box* b, int L );
    
    // Secondary algorithm to determine transfer pairs and close pairs
    void defineTransferAndClose( Box* b1, Box* b2, int L );
    
  };
  
}

#endif
