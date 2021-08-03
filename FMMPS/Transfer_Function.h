/*
 *  Transfer_Function.h
 *  PSCAcoustic
 *
 *  Transfer between boxes at given level (high- and low-frequency)
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

#ifndef TRANSFER_FUNCTION_FMM_H
#define TRANSFER_FUNCTION_FMM_H

#include "Quadrature.h"
#include "TransferUtils.h"

namespace FMM {
  

  //---- Super class for transfer function ----//
  class Transfer_Function
  {
    
  public:

    static Indexing global_index;

    // Constructor/Destructor
    Transfer_Function(){}
    virtual ~Transfer_Function(){}
    
    // Execution
    // Transfer from vec1 to vec2
    virtual void Transfer( std::vector<complex>& vec1, std::vector<complex>& vec2, Type type ) { return; }  
  };
  
  
  
  //---- High-frequency transfer function ----//
  class HF_Transfer_Function : public Transfer_Function
  {
    std::vector<complex> C;
    complex k_out;

  public:
    
    // Constructor
    HF_Transfer_Function(Quadrature* quad, const Vec3& r0, complex k_out_);
    
    //Destructor
    ~HF_Transfer_Function(){}

    // Execution
    // from vecFrom to vecTo
    virtual void Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type );
  };




  //---- Low-frequency transfer function ----//
  class Direct_Transfer_Function : public Transfer_Function 
  {
 
    PointAndShoot* PS_Transfer;

  public:
    // Constructor/Destructor
    //Direct_Transfer_Function( Indexing* index, Vec3 r, complex k_out);
    Direct_Transfer_Function( Indexing* index, Vec3 r, complex k_out);
    ~Direct_Transfer_Function();

    // From vecFrom to vecTo
    virtual void Transfer( std::vector<complex>& vecTo, std::vector<complex>& vecFrom, Type type );
  };

}

#endif
