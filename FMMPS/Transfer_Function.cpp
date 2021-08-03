/*
 *  Transfer_Function.cpp
 *  PSCAcoustic
 *
 *  Transfer between boxes at given level (high- and low-frequency)
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

#ifndef TRANSFER_FUNCTION_FMM_CPP
#define TRANSFER_FUNCTION_FMM_CPP

#include "General.h"
#include "../Coordinates.h"
#include "HelmholtzUtil.h"
#include "TransferUtils.h"
#include "Transfer_Function.h"

namespace FMM {

  Indexing Transfer_Function::global_index(LMAX);

  //---- High-frequency transfer function ----//
  
  HF_Transfer_Function::HF_Transfer_Function( Quadrature* quad, const Vec3& r0, complex k_out_) 
    : Transfer_Function(), k_out(k_out_)
  {
    int L = quad->L;
    int N_rows = quad->numRows();
    int N_phi = 2*L + 2;
    C.resize( quad->size() );
    
    
    // Get the low-pass modified transfer matrix N_rows x N_phi 
    std::vector< std::vector<complex> > T = mod_transfer_fn(L, k_out, N_rows, N_phi, r0);
    
    // maxRowSize can be larger than Nphi -> shouldn't interpolate in place
    std::vector<complex> t( std::max(quad->maxRowSize(),N_phi) );
    
    // Anterpolate the rows for an optimized quadrature
    for( int n = 0; n < N_rows; ++n ) {
      QuadratureRow& row = quad->getRow(n);
      int N_phi_n = row.size();
      
      // Anterpolate to N_phi_k
      for( int m = 0; m < N_phi; ++m ) t[m] = T[n][m];
      fftinterp( &t[0], N_phi, N_phi_n );
      
      // Copy into the NFunction
      for( int m = 0; m < N_phi_n; ++m ) {
	int index = row.getPoint(m).index;
	C[index] = t[m];
      }
    }
  }  
  
  
  
  // Execution
  // from vecFrom to vecTo
  void HF_Transfer_Function::Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type ) 
  { 
    assert( vecTo.size() == C.size() );
    assert( vecTo.size() == vecFrom.size() );
    
    int N = vecFrom.size();
    // ***This factor comes from the definition of Cecka & Darve which uses a factor CI*K_OUT/(4*PI) in 
    // the definition of the the transfer function. My implementation does not use such factor. It must
    // therefore be removed. Maybe it could be adressed directly in the definition of the transfer function...
    //if( type == FORWARD ){
      for( int i = 0; i < N; i++ )
	vecTo[i] += C[i] * vecFrom[i] * 4*PI/(k_out * CI);
      //} else if( type == ADJOINT ){
      //for( int i = 0; i < N; i++ )
      //	vecTo[i] += conj(C[i]) * vecFrom[i] * conj(4*PI/(k_out * CI)); //TODO : Is normalization ok?
      //}
    
  }
  
  
  
  

  //---- Low-frequency transfer function ----//
  
  Direct_Transfer_Function::Direct_Transfer_Function( Indexing* index, Vec3 r, complex k_out) : Transfer_Function()    
  {
    double rad, theta, phi;
    Cart2Sf(r.x, r.y, r.z, rad, theta, phi);
    
    // Initialize P&S transfer
    PS_Transfer = new PointAndShoot(rad, k_out, theta, phi, index, index, false);
  }
  
  //Destructor
  Direct_Transfer_Function::~Direct_Transfer_Function()
  { 
    delete PS_Transfer; 
  }
  
  // Execution
  // From vecFrom to vecTo
  void Direct_Transfer_Function::Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type ) 
  { 
    assert( (int) vecFrom.size() == (int) vecTo.size() );
    std::vector<complex> vec = PS_Transfer->Apply( vecFrom, type );
    assert( (int) vec.size() == (int) vecTo.size() );
    
    // Update expansion
    for( int k = 0; k < (int) vecTo.size(); k++ )
      vecTo[k] += vec[k];
  }


}

#endif
