/*
 *  Transfer_Function.h
 *  PSCAcoustic
 *
 *  Transfer between boxes at given level (high- and low-frequency)
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
