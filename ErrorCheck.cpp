#ifndef ERRORCHECK_CPP
#define ERRORCHECK_CPP

/*
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved. 
 *
 */


#include <vector>
#include <complex>
#include <iostream>

#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"
#include "Distributions.h"
#include "PSCEnv.h"
#include "Imaging.h"
#include "Solver.h"
#include "Finalize.h"

#include "./FMMPS/Vec3.h"
#include "./FMMPS/Quadrature.h"
#include "./FMMPS/TransferUtils.h"
#include "./FMMPS/Transfer_Function.h"
#include "./FMMPS/Translation_Function.h"
#include "./FMMPS/SphericalHarmonicsTransform.h"
#include "./L_choice.h"
#include "./FMMPS/HelmholtzUtil.h"


int main(int argc,char **args)
{



// Global indexing
  Indexing Index(LMAX); // TODO: Should not have different instances of Index
  FMM::HelmKernel KERNEL;
  KERNEL.kappa = K_OUT;
  
  double EPS;
  bool close = true;
  bool LF = true;   // low-frequency or high-frequency check
  int N = 100;
  double boxSize = 1.;
  double R = 1.1*RADIUS;
  double d_min = 0.1 * 2*PI/std::abs(K_OUT);
  std::string filename("../AcousticOutputFiles/ErrorCheck_boxSize=1_LF_A.csv");
  
 // *** Scatterer clusters *** 
 Pvec center_1(0.,0.,0.,Pvec::Cartesian);
 Pvec center_2(sqrt(3.)*boxSize, sqrt(3.)*boxSize, sqrt(3.)*boxSize,Pvec::Cartesian);
 //Pvec center_2(3*boxSize, 0., 0.,Pvec::Cartesian);

 cout << "***Building scatterers...";
 std::vector<Scatterer> ScL_1;
 std::vector<Pvec> scat_loc = RandSpherical(boxSize, RADIUS, D_MIN, center_1, N);
 assert( N == (int) scat_loc.size() );
 for( int i = 0; i < N; i++ ){
   Scatterer scatterer(RADIUS, K, K_OUT, RHO, scat_loc[i]);
   ScL_1.push_back(scatterer);
   //cout << scat_loc[i].x << " : " << scat_loc[i].y << " : " << scat_loc[i].z << " : " << endl;
 }
 
 std::vector<Scatterer> ScL_2;
 scat_loc = RandSpherical(boxSize, RADIUS, D_MIN, center_2, N);
 assert( N == scat_loc.size() );
 for( int i = 0; i < N; i++ ){
   Scatterer scatterer(RADIUS, K, K_OUT, RHO, scat_loc[i]);
   ScL_2.push_back(scatterer);
   //cout << scat_loc[i].x << " : " << scat_loc[i].y << " : " << scat_loc[i].z << " : " << endl;
 }
 cout << "Done" << endl;

 



 cout << "***Computing right-hand side...";
 Pvec dir(0.,0.,1.,Pvec::Cartesian);
 PlaneWave IW(1.+0.*CI, dir);
 std::vector< std::vector<complex> > RHS( N, std::vector<complex>((LMAX+1)*(LMAX+1)) );
 std::vector<complex> rhs((LMAX+1)*(LMAX+1));
 // For each scatterer in box 1 compute r.h.s.
 for( int n = 0 ; n < N ; n++ )
   {
     //Transfer expansion about center of scatterer
     Pvec vec = ScL_1[n].getLocation();
     IW.Translate(vec, rhs);
     
     //Apply T-matrix
     ScL_1[n].TM_Apply(rhs);
     
     for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ )
       RHS[n][i] = rhs[i];

     
   }
 cout << " Done " << endl;
 
 
   // Initialiaze static variables for Transfer Utils
   FMM::S_Rotation::S_Init( std::ceil(2.*LMAX) );
   FMM::Z_transfer::Z_Init( std::ceil(2.*LMAX) );


   cout << "***Computing exact solution...";
   std::vector< std::vector<complex> > EXP_1( N, std::vector<complex>((LMAX+1)*(LMAX+1)) );
   std::vector< std::vector<complex> > EXP_2( N, std::vector<complex>((LMAX+1)*(LMAX+1)) );
   // Perform transfer
   std::vector<complex> v;
   for( int i = 0; i < N; i++ ){
     std::vector<complex> TM = ScL_2[i].getTMatrix();
     for( int j = 0; j < N; j++ ){
       Pvec vec = ScL_2[i].getLocation() - ScL_1[j].getLocation();
       
       FMM::PointAndShoot transfer(vec.r, K_OUT, vec.theta, vec.phi, &Index, &Index, false);
       v = transfer.Apply( RHS[j], FMM::FORWARD );
       for( int k = 0; k < Index.size(); k++){
	 int l = Index(k,0);
	 EXP_2[i][k] += TM[l] * v[k];
       }
       
       if( (i != j) && close){
	 Pvec vec = ScL_1[i].getLocation() - ScL_1[j].getLocation();
	 FMM::PointAndShoot transfer_1(vec.r, K_OUT, vec.theta, vec.phi, &Index, &Index, false);
	 v = transfer_1.Apply( RHS[j], FMM::FORWARD );
	 
	 for( int k = 0; k < Index.size(); k++){
	   int l = Index(k,0);
	   EXP_1[i][k] += TM[l] * v[k];
	 }
       }
       
     }
   }
   cout << " Done " << endl;
   


 int M = 9;
 std::vector<int> L_S(M);
 std::vector<int> L_T(M);
 std::vector<int> Q_SIZE(M);
 std::vector<double> Error_S(M);
 std::vector<double> RelError_S(M);
 std::vector<double> Error_T(M);
 std::vector<double> RelError_T(M);
 std::vector<double> eps(M);
 int maxIndexSize = 2.*LMAX;
 for( int p = 0; p < M; p ++ ){
   FMM::Quadrature* quad;
   EPS = pow(10., -3. - (double) p );
   eps[p] = EPS;
   
   cout << "Computing parameters..." << endl;
   L_S[p] = L_max_scat(boxSize, d_min, eps[p], K, K_OUT);
   cout << "L_S : " << L_S[p] << endl;
   L_T[p] = L_max_level(boxSize, eps[p], RADIUS, K, K_OUT);
   cout << "L_T : " << L_T[p] << endl;

   if( !LF )
     quad = new FMM::Quadrature(KERNEL, boxSize, EPS, K, K_OUT);
   

   cout << "Done: Computing parameters" << endl << endl;
   
   cout << "EPS : " << EPS << " ; L_S : " << L_S[p] << endl;
   cout << "EPS : " << EPS << " ; L_T : " << L_T[p] << endl << endl;
   
   Indexing Index_S(L_S[p]);
   Indexing Index_T(L_T[p]);
   
   // Initialiaze static variables for Transfer Utils
   int tempIndexSize = LMAX;
   tempIndexSize = std::max(tempIndexSize, L_S[p]);
   tempIndexSize = std::max(tempIndexSize, L_T[p]);

   // Initialiaze static variables for Transfer Utils
   if( tempIndexSize > maxIndexSize ){
     FMM::S_Rotation::S_Delete();
     FMM::S_Rotation::S_Init( std::ceil(1.2*tempIndexSize) );
     FMM::Z_transfer::Z_Delete();
     FMM::Z_transfer::Z_Init( std::ceil(1.2*tempIndexSize) );
     maxIndexSize = tempIndexSize;
   }
   
     //cout << "max size : " << maxIndexSize << endl;
   if( !LF )
     Q_SIZE[p] = quad->size();
   
   
   
   
   

   
   // ***!!! Transfer from box 1 to box 2 and check error
   
   
    // Right-hand side
   if( close ){
     std::vector< std::vector<complex> > RHS_S( N, std::vector<complex>((L_S[p]+1)*(L_S[p]+1)) );
     std::vector<complex> rhs((LMAX+1)*(LMAX+1));
     
     // For each scatterer in box 1 compute r.h.s.
     for( int n = 0 ; n < N ; n++ )
       for( int i = 0; i < Index_S.size(); i++ )
	 RHS_S[n][i] = RHS[n][i];
     
     
     // Perform transfer
     std::vector< std::vector<complex> > EXP_S( N, std::vector<complex>((L_S[p]+1)*(L_S[p]+1)) );
     for( int i = 0; i < N; i++ ){
       std::vector<complex> TM = ScL_2[i].getTMatrix();
       for( int j = 0; j < N; j++ ){
	 
	 if( i != j ){
	   Pvec vec = ScL_1[i].getLocation() - ScL_1[j].getLocation();
	   
	   FMM::PointAndShoot transfer_2(vec.r, K_OUT, vec.theta, vec.phi, &Index_S, &Index_S, false);
	   std::vector<complex> v = transfer_2.Apply( RHS_S[j], FMM::FORWARD );
	   
	   for( int k = 0; k < Index_S.size(); k++){
	     int l = Index_S(k,0);
	     EXP_S[i][k] += TM[l] * v[k];
	   }
	 }
	 
       }
     }
     
     
     // Compute error
     Error_S[p] = 0;
     double Norm_S = 0;
     for( int i = 0; i < N; i++ ){
       for( int k = 0; k < Index_S.size(); k++){
	 int l = Index_S(k,0);
	 
	 Error_S[p] += std::abs( (EXP_1[i][k] - EXP_S[i][k]) *  gsl_sf_hankel_1(l, std::real(K_OUT)*R) );
	 Norm_S += std::abs( EXP_1[i][k] *  gsl_sf_hankel_1(l, std::real(K_OUT)*R) );
       }
     }
     RelError_S[p] = Error_S[p] / Norm_S;
     
     cout << "EPS : " << EPS << endl;
     cout << "absolute Error_S : " << Error_S[p] << endl;
     cout << "relative Error_S : " << RelError_S[p] << endl << endl;
     
   }
   

   
   
   
   
   // *** Testing L_T for scatterers only *** 

   // Perform transfer
   std::vector< std::vector<complex> > EXP_T(N, std::vector<complex>((LMAX+1)*(LMAX+1)) );
   std::vector<complex> EXP_T_1;
   std::vector<complex> EXP_T_2;
   
   if( LF ){
     EXP_T_1.resize( (L_T[p]+1)*(L_T[p]+1) );
     EXP_T_2.resize( (L_T[p]+1)*(L_T[p]+1) );
   } else {
     EXP_T_1.resize( quad->size() );
     EXP_T_2.resize( quad->size() );
     for( int q = 0; q < quad->size(); q++){
       EXP_T_1[q] = 0.;
       EXP_T_2[q] = 0.;
     }
   }
 

   for( int i = 0; i < N; i++ ){
     Pvec vec_1 = center_1 - ScL_1[i].getLocation();

     if( LF ){
       std::vector<complex> v;
       FMM::PointAndShoot transfer_1(vec_1.r, K_OUT, vec_1.theta, vec_1.phi, &Index_T, &Index, true);
       v = transfer_1.Apply( RHS[i], FMM::FORWARD );

       for( int k = 0; k < Index_T.size(); k++)
	 EXP_T_1[k] += v[k];
     } else {
       std::vector<complex> v(quad->size());
       for( int q = 0; q < quad->size(); q++ )
	 v[q] = 0.;
       FMM::Vec3 v3(vec_1.x, vec_1.y, vec_1.z);
       FMM::Leaf_HF_Translation_Function transfer_1(v3, K_OUT, quad);
       transfer_1.TranslateUp(v, RHS[i], FMM::FORWARD); 
       
       for( int k = 0; k < quad->size(); k++)
	 EXP_T_1[k] += v[k];
     }

   }

   Pvec vec_t = center_2 - center_1;
   if( LF ){
     FMM::PointAndShoot TransferCluster(vec_t.r, K_OUT, vec_t.theta, vec_t.phi, &Index_T, &Index_T, false);
     EXP_T_2 = TransferCluster.Apply( EXP_T_1, FMM::FORWARD );
   } else {
     FMM::Vec3 v3(vec_t.x, vec_t.y, vec_t.z);
     FMM::HF_Transfer_Function TransferCluster(quad, v3, K_OUT);
     TransferCluster.Transfer( EXP_T_2, EXP_T_1, FMM::FORWARD);
   }
   
   for( int j = 0; j < N; j++ ){
     Pvec vec_2 = ScL_2[j].getLocation() - center_2;
     std::vector<complex> w((LMAX+1)*(LMAX+1));

     if( LF ){
       FMM::PointAndShoot transfer_2(vec_2.r, K_OUT, vec_2.theta, vec_2.phi, &Index, &Index_T, true);
       w = transfer_2.Apply( EXP_T_2, FMM::FORWARD );
     } else {
       for( int q = 0; q < w.size(); q++ )
	 w[q] = 0.;
       FMM::Vec3 v3(-vec_2.x, -vec_2.y, -vec_2.z);
       FMM::Leaf_HF_Translation_Function transfer_2(v3, K_OUT, quad);
       transfer_2.TranslateDown(w, EXP_T_2, FMM::FORWARD); 
     }
     
     // Apply T-matrix
     ScL_2[j].TM_Apply(w);
     
     for( int k = 0; k < Index.size(); k++)
       EXP_T[j][k] = w[k];
   }
   
   
   // Compute error
   Error_T[p] = 0;
   double Norm_T = 0;
   for( int i = 0; i < N; i++ ){
     for( int k = 0; k < (LMAX+1)*(LMAX+1); k++){
       int l = Index(k,0);
       Error_T[p] += std::abs( (EXP_2[i][k] - EXP_T[i][k]) *  gsl_sf_hankel_1(l, std::real(K_OUT)*R) );
       Norm_T += std::abs( EXP_2[i][k] *  gsl_sf_hankel_1(l, std::real(K_OUT)*R) );

       //cout << EXP_2[i][k] << " : " << EXP_T[i][k] << endl;
     }
   }
   RelError_T[p] = Error_T[p] / Norm_T;
   
   cout << "EPS : " << EPS << endl;
   cout << "absolute Error_T : " << Error_T[p] << endl;
   cout << "relative Error_T : " << RelError_T[p] << endl << endl;
   
 }


 ofstream file((char*) filename.c_str(), ios::out);
 for( int p = 0; p < M; p++ ){
   if( LF ){
     file << std::setprecision(15) << eps[p] << "," << L_S[p] << "," << Error_S[p] << "," << RelError_S[p] << "," <<
       L_T[p] << "," <<  Error_T[p] << "," << RelError_T[p] << endl;
   } else {
     file << std::setprecision(15) << eps[p] << "," << L_S[p] << "," << Error_S[p] << "," << RelError_S[p] << "," <<
       Q_SIZE[p] << "," <<  Error_T[p] << "," << RelError_T[p] << endl;
   }
 }
 file.close();


  return 0;
}





#endif
