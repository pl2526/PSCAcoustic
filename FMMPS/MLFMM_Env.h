/*
 *  MLFMM_Env.h
 *  PSCAcoustic
 *
 *  Multi-level FMM environment
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
