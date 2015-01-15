#ifndef LEVEL_FMM_CPP
#define LEVEL_FMM_CPP

/*
 *  Level.cpp
 *  PSCAcoustic
 *
 *  Objects needed at each level of the FMM tree
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

#include "Level.h"

namespace FMM {

  Indexing Level::global_index(LMAX);

  // Constructors
  Level::Level( double boxSize_, complex k_out_, complex k_, double cutoff) 
    : boxSize(boxSize_), k_out(k_out_), k(k_),
      quad(NULL), translationFunc(BRANCH),
      upReterp(NULL), downReterp(NULL),
      SH_T(NULL)
  {
    

    // TODO: SHOULD BE DIFFERENT FOR EACH
    //double c = std::min(std::abs(K_L_OUT), std::abs(K_T_OUT));
    //HF = ( boxSize > std::abs(CUTOFF_SIZE/(2*PI*OMEGA/c)) ) ? true : false;

    HF = ( boxSize > CUTOFF_SIZE) ? true : false;
    INTF = false;
    LEAF = false;
    
    if( verbose >= 2){
      cout << "boxSize : " << boxSize_ << endl;
      cout << "HF : " << HF << endl;
    }
    
  }
  
  
  
  // Destructor
  Level::~Level() 
  {
    
    // Delete transfers
    // Only delete those assigned once. Deleting copies creates double free errors.
    Trans_Vector LastTV;
    for( TVIter tvi = transList.begin(); tvi != transList.end(); ++tvi ) {
      Trans_Vector& tv = *tvi;
      
      if( LastTV.x != tv.x || LastTV.y != tv.y || LastTV.z != tv.z )
	delete tv.T;   
      
      LastTV = tv;
    }
    
    
    // Delete Translation Functions
    for( TLIter tli = translationFunc.begin();
	 tli != translationFunc.end();
	 ++tli ) {
      Translation_Function* TL = *tli;
      delete TL;
    }
    
    
    // Delete close transfers
    for( TFIter ctf = closeFunc.begin(); 
	 ctf != closeFunc.end();
	 ++ctf ){
      Transfer_Function* T = *ctf;
      delete T;
    }
    
    
    if( SH_T != NULL )
      delete SH_T;
    delete quad;
    delete upReterp;
    delete downReterp;
  }
  
  
  
  //-- Auxilliary functions --//      
  
  void Level::setL(double eps, double r, complex k, complex k_out)
  {


    // TODO: PASS RIGHT K
    L =  L_max_level(boxSize, eps, r, K, k_out);

    
    if( !HF )
      index = new Indexing(L);

    cout << " level L : " << L << endl;
  }
  
  void Level::initializeFields( Type type )
  {
    // Zero the box multipole and locals
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ) {
      Box* b = *bi;
      
      if( HF ) {
	b->makeMultipole( quad->size() );
	b->makeLocal( quad->size() );
      } else {
	b->makeMultipole( index->size() );
	b->makeLocal( index->size() );
      }
      
    }
  }
  
  void Level::zeroFields()
  {
    // Zero the box multipole and locals
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ) {
      Box* b = *bi;
      b->getMultipole().assign( b->getMultipole().size(), 0.);
      b->getLocal().assign( b->getLocal().size(), 0.);
    }
  }
  
  
  
  
  //--- Quadrature ---//
  
  // TODO : Should spherical harmonics transform be constructed at the same time as quad?
  void Level::defineQuadrature( HelmKernel KERNEL, double eps)
  {
    // Construct the Quadrature
    quad = new Quadrature(KERNEL, boxSize, eps, k, k_out);
  }
  
  

  //--- Spherical Harmonics Transform ---//
  
  void Level::computeSH_T_Forward( Type type )
  {
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ){
      Box* b = *bi;
      
      if( type == FORWARD )
	b->getMultipole() = SH_T->Forward( b->getMultipole() );
      else if( type == ADJOINT )
	b->getLocal() = SH_T->Forward( b->getLocal() );
      
    }
  }
  
  void Level::computeSH_T_Backward( Type type )
  {
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ){
      Box* b = *bi;
      
      if( type == FORWARD ){
	b->getLocal() = SH_T->Backward( b->getLocal() );
      } else if( type == ADJOINT ){
	b->getMultipole() = SH_T->Backward( b->getMultipole() );
      }
      
      
    }
  }
  
  
  
  
  //--- Translation and Transfer ---//
  
  void Level::defineTransfers( double eps, complex k_out )
  {
    transList.sort();
    
    int k = 0;	
    
    Trans_Vector LastTV;
    for( TVIter tvi = transList.begin(); tvi != transList.end(); ++tvi ) {
      Trans_Vector& tv = *tvi;
      
      if( LastTV.x != tv.x || LastTV.y != tv.y || LastTV.z != tv.z ) {
	// Compute a new Transfer Function
	Vec3 r0( tv.x * boxSize, tv.y * boxSize, tv.z * boxSize );
	
	if( HF ) {
	  tv.T = new HF_Transfer_Function( quad, r0, k_out );
	} else {
	  tv.T = new Direct_Transfer_Function( index, r0, k_out );
	}
	k++;
      } else {
	//This Transfer Function is the same as the last one
	tv.T = LastTV.T;
      }
      
      LastTV = tv;
    }
    
    cout << "Number of transfers defined : " << k << endl;
  }
  

  void Level::addTransfer( Box* b1, Box* b2 )
  {
    Vec3 r0 = b2->center - b1->center;
    r0 /= boxSize;
    r0.x = round(r0.x); r0.y = round(r0.y); r0.z = round(r0.z);
    
    transList.push_back( Trans_Vector(b1, b2,  r0) );
    transList.push_back( Trans_Vector(b2, b1, -r0) );
  }
 
  
  void Level::addClose( Box* b1, Box* b2 )
  { 
    closeList.push_back( Close_Vector(b1, b2) );
  }
  
  
  // Define close transfers (only at leaf level)
  void Level::addCloseTransfer( int& idxTo, int& idxFrom, Vec3& r0, complex k_out )
  {
    std::vector<int> idxVec(2);
    idxVec[0] = idxTo;
    idxVec[1] = idxFrom;
    closeIdx.push_back( idxVec );
    closeVec.push_back(r0);
    
    if( CT_PRECOMPUTE ){
      Direct_Transfer_Function* closeTransfer = new Direct_Transfer_Function( &global_index, r0, k_out);
      closeFunc.push_back( closeTransfer );
    }
    
  }
  
  
  // Define the translation functions from a lower level
  void Level::defineTranslations(Indexing* pIndex, Indexing* cIndex, complex k_out)
  {
    // Depending on BRANCH, construct the translation std::vectors
    for( int k = 0; k < BRANCH; ++k ) {
      // Generate the kth translation std::vector
      // Could do this alot better with the NTree...
      Vec3 r;
      if( DIM == 1 ) {
	if( k == 0 ) r = Vec3( boxSize/4, 0, 0);
	if( k == 1 ) r = Vec3(-boxSize/4, 0, 0);
      }
      if( DIM == 2 ) {
	if( k == 0 ) r = Vec3( boxSize/4, boxSize/4, 0);
	if( k == 1 ) r = Vec3( boxSize/4,-boxSize/4, 0);
	if( k == 2 ) r = Vec3(-boxSize/4, boxSize/4, 0);
	if( k == 3 ) r = Vec3(-boxSize/4,-boxSize/4, 0);
      }
      if( DIM == 3 ) {
	if( k == 0 ) r = Vec3( boxSize/4, boxSize/4, boxSize/4);
	if( k == 1 ) r = Vec3( boxSize/4, boxSize/4,-boxSize/4);
	if( k == 2 ) r = Vec3( boxSize/4,-boxSize/4, boxSize/4);
	if( k == 3 ) r = Vec3( boxSize/4,-boxSize/4,-boxSize/4);
	if( k == 4 ) r = Vec3(-boxSize/4, boxSize/4, boxSize/4);
	if( k == 5 ) r = Vec3(-boxSize/4, boxSize/4,-boxSize/4);
	if( k == 6 ) r = Vec3(-boxSize/4,-boxSize/4, boxSize/4);
	if( k == 7 ) r = Vec3(-boxSize/4,-boxSize/4,-boxSize/4);
      }
      
      
      // Construct translation based on frequency regime
      if( HF )
	translationFunc[k] = new HF_Translation_Function(r, k_out, INTF, SH_T, quad);
      else 
	translationFunc[k] = new Direct_Translation_Function(r, k_out, pIndex, index, cIndex);
      
    }
  }
  
  
  
  //--- Execution ---//
  
  void Level::applyTransferFunctions( Type type )
  {
    for( TVIter tvi = transList.begin(); tvi != transList.end(); ++tvi ) {     
      
      Transfer_Function* T = tvi->T;
      
	  Box* b1 = tvi->b1;
	  Box* b2 = tvi->b2;
	  
	  if( type == FORWARD ){
	    // Transfer from box 1(multipole) to box 2(local)
	    T->Transfer( b2->getLocal(), b1->getMultipole(), FORWARD);
	  } else if( type == ADJOINT ){
	    // Adjoint transfer from box 2(local) to box 1(multipole)
	    T->Transfer( b1->getMultipole(), b2->getLocal(), ADJOINT);
	  }
	  
    }
  }
  
  
  void Level::applyClose(const std::vector< std::vector<complex> >& PWE, 
			 std::vector< std::vector<complex> >& T_PWE, 
			 complex k_out, Type type )
  {
    int N = closeIdx.size();
    for( int n = 0; n < N; n++ ){
      std::vector<complex> vecTo(global_index.size());
      int idxTo = 0;
      int  idxFrom = 0;
      
      if( type == FORWARD ){
	idxTo = closeIdx[n][0];
	idxFrom = closeIdx[n][1];
      }else if( type == ADJOINT ){
	idxTo = closeIdx[n][1];
	idxFrom = closeIdx[n][0];
      }
      
      
      // Proceed to transfer
      std::vector<complex> vecFrom = PWE[idxFrom]; 
      if( CT_PRECOMPUTE ){
	closeFunc[n]->Transfer(vecTo, vecFrom, type);
      } else {
	Direct_Transfer_Function* closeTransfer = new Direct_Transfer_Function( index, closeVec[n], k_out);
	closeTransfer->Transfer(vecTo, vecFrom, type);
	delete closeTransfer;
      }
      
      // Update expansion
      for( int i = 0; i< global_index.size(); i++ )
	T_PWE[idxTo][i] += vecTo[i];
      
    }
    
  }
  
  
}
#endif
