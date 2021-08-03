/*
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

#ifndef TU_FMM_CPP
#define TU_FMM_CPP

#include "TransferUtils.h"

namespace FMM {


  //**** Point and Shoot Implementation ****//
  
  //---Rotation about z-axis---//

  // Constructor
  Z_Rotation::Z_Rotation( Indexing* index, double theta ){ Compute(matrix, index, theta); }
  
  void Z_Rotation::Compute( std::vector<complex>& matrix, Indexing* index, double theta ){
    
    matrix.resize( index->size() );    
    
    int m;
    complex entry;
    for( int i = 0; i < (int) index->size(); i++ ){
      m = (*index)(i,1);
      entry = exp( m * theta * CI );
      
      matrix[i] = entry;
    }
    
  }
  
  // Proceed to translation in-place
  void Z_Rotation::Apply( std::vector<complex>& vec, Indexing* index, Sign S ){
    
    int N = matrix.size();
    int M = index->size();

    assert( (int) vec.size() == M );
    assert( M <= N );
    
    if( S == POS ){
      for( int n = 0; n < M; n++)
	vec[n] *= matrix[n];      
    } else {
      for( int n = 0; n < M; n++)
	vec[n] *= conj( matrix[n] );     
    }
    
  }
  
  
  
  
  
  
  
  //--- S Matrix ---//

    
    
    
  void S_Rotation::S_Init( int L ) {
    S_Rotation::Init = true;
    // Only nonzero if l's are the same
    // from (l,m2) to (l,m1)
    index = new Indexing(L);
    
    //TODO : this matrix embeds theta=pi/2 only
    complex D[3][3];
    D[0][0] = 1./2.;
    D[0][1] = 1./sqrt(2);
    D[0][2] = 1./2.;
    D[1][0] = -1./sqrt(2);
    D[1][1] = 0.;
    D[1][2] = 1./sqrt(2);
    D[2][0] = 1./2.;
    D[2][1] = -1./sqrt(2);
    D[2][2] = 1./2.; 
    
    
    //Allocate space for matrix storage (in sparse form)
    S_matrix = new complex*[(L+1)*(L+1)];
    for( int k = 0; k < (L+1)*(L+1); k++ ){
      S_matrix[k] = new complex[(L+1)*(L+1)];
    }
    
    // Populate matrix entries through recursion (based on Choi et al. paper)
    S_matrix[0][0] = 1.;
    for( int m_p = 0; m_p <= L; m_p++ ){
      
      // Case  l = m_p
      complex *s;
      if( m_p != 0 ){
	int l = m_p;
	for( int m = -l; m <= l; m++ ){
	  s = &S_matrix[index->LocateIndex(l,m_p)][index->LocateIndex(l,m)];
	  *s = 0;
	  
	  if( (l-1) >= std::abs(m) )
	    *s += c(l,m,m_p) * D[1][2] * S_matrix[index->LocateIndex(l-1,m_p-1)][index->LocateIndex(l-1,m)];
	  
	  if( (l-1) >= std::abs(m-1) )
	    *s += d(l,m,m_p) * D[2][2] * S_matrix[index->LocateIndex(l-1,m_p-1)][index->LocateIndex(l-1,m-1)];
	  
	  if( (l-1) >= std::abs(m+1) )
	    *s += d(l,-m,m_p) * D[0][2] * S_matrix[index->LocateIndex(l-1,m_p-1)][index->LocateIndex(l-1,m+1)];
	  
	  // Negative index through symmetry
	  S_matrix[index->LocateIndex(l,-m_p)][index->LocateIndex(l,-m)] = pow(-1., (double)(m-m_p)) * (*s);
	  
	}
      }
      
      
      // Case l > m_p
      for( int l = (m_p+1); l <= L; l++ ){
	for( int m = -l; m <=l; m++ ){
	  
	  s = &S_matrix[index->LocateIndex(l,m_p)][index->LocateIndex(l,m)];
	  *s = 0;
	  
	  if( (l-1) >= std::abs(m) && (l-1) >= std::abs(m_p) )
	    *s += a(l,m,m_p) * D[1][1] * S_matrix[index->LocateIndex(l-1,m_p)][index->LocateIndex(l-1,m)];
	  
	  if( (l-1) >= std::abs(m-1) && (l-1) >= std::abs(m_p) )
	    *s += b(l,m,m_p) * D[2][1] * S_matrix[index->LocateIndex(l-1,m_p)][index->LocateIndex(l-1,m-1)];
	  
	  if( (l-1) >= std::abs(m+1) && (l-1) >= std::abs(m_p) )
	    *s += b(l,-m,m_p) * D[0][1] * S_matrix[index->LocateIndex(l-1,m_p)][index->LocateIndex(l-1,m+1)];
	  
	  // Negative index through symmetry
	  S_matrix[index->LocateIndex(l,-m_p)][index->LocateIndex(l,-m)] = pow(-1., (double)(m-m_p)) * (*s);
	}
      }
      
    }
    
    
    // TODO : should be stored in sparse form
    // Compute entry for D(0,Pi/2,0) * D(Pi/2,0,0)
    for( int l_p = 0; l_p <= L; l_p++ ){	
      for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	for( int l = 0; l <= L; l++ ){
	  for( int m = -l; m <= l; m++ ){
	    
	    S_matrix[index->LocateIndex(l_p,m_p)][index->LocateIndex(l,m)] *= exp( m * PI/2 * CI);
	    
	  }
	}
      }
    }
    
    
  }
  


  // Delete static pointer
  void S_Rotation::S_Delete(){
    int K = index->size();
    delete index;
    
    for( int k = 0; k < K; k++ )
      delete S_matrix[k]; 
    delete S_matrix;
  }
    
    
    
    
  // Apply Rotation operator
  void S_Rotation::Apply( std::vector<complex>& vec, int L, Type type ){
    
    int idx1, idx2;
    std::vector<complex> V( (L+1)*(L+1) );
    for( int l = 0; l <= L; l++ ){	
      for( int m1 = -l; m1 <= l; m1++ ){
	for( int m2 = -l; m2 <= l; m2++ ){
	  
	  idx1 = index->LocateIndex(l, m1);
	  idx2 = index->LocateIndex(l, m2);
	  
	  if( type == ADJOINT ){
	    V[idx1] += conj(S_matrix[idx2][idx1]) * vec[idx2];
	  } else if( type == FORWARD ){
	    V[idx1] += S_matrix[idx1][idx2] * vec[idx2];
	  }
	  
	}	
      }
    }
    
    
    // TODO : Do I need this loop? Will compiler unwrap it?
    for( int k = 0; k < (int) vec.size(); k++ )
      vec[k] = V[k];
    
  }

  complex** S_Rotation::S_matrix;
  bool S_Rotation::Init = false;
  Indexing* S_Rotation::index;
  
  





    //--- Transfer along z-axis ---//
    

  Z_transfer::Z_transfer( double R, complex k_out, Indexing* index1_, Indexing* index2_, bool reg)
      : index1(index1_), index2(index2_)
  {
    //cout << "R : " << R << endl;
    assert( R > 0 );
    
    // Compute all transfer coefficients for degree <= L
    // (only for transfer/irregular)
    
    // Identify index corresponding to nonzero and store in sparse form
    // (Only nonzero transfer coefficient if m's are the same)
    // vec_idx[0] : l1
    // vec_idx[1] : m1
    // vec_idx[2] : l2
    // vec_idx[3] : m2
    std::vector<int> vec_idx(4);
    complex entry;
    
    // TODO : clean up L_1 and L_2
    L_1 = (*index1_)(index1_->size()-1, 0);
    L_2 =  (*index2_)(index2_->size()-1, 0);
    int L_max = std::max(L_1, L_2);
    
    // Compute coefficients
    // (use sparse matrix for construction but store in std::vector for speed and memory)
    SpMatrix coeff2;
    Compute(coeff2, L_max, R, k_out, reg);
    
    L = std::min(L_1, L_2);
    for( int m = -L; m <= L; m++ )
      for( int l1 = std::abs(m); l1 <= L_1; l1++ )
	for( int l2 = std::abs(m); l2 <= L_2; l2++ )
	  coeff.push_back( coeff2(index1->LocateIndex(l1,m), index2->LocateIndex(l2,m)) );
    
  }
  
  
  // Initialization for translation along z-axis
  void Z_transfer::Z_Init( int L ){ 
    Z_transfer::Init = true;
    index = new Indexing(L); 
  }

  // Destruction of translation along z-axis
  void Z_transfer::Z_Delete(){ delete index; }

      
    // Recursive computation of transfer  and translation coefficients
    // For details, see Gumerov, Duraiswami : "Fast, Exact, and Stable Computation of Multipole Translation
    // and Rotation Coefficients for the 3-D Helmholtz equation"
    // Matlab Prototype : GumerovTransferCoefficients.m
  double a(int n,int m){
    double entry;
    ( n >= std::abs(m) ) ? entry = sqrt( (n+1.+std::abs(m)) * (n+1.-std::abs(m)) / ((2.*n+1.)*(2.*n+3.)) ) : entry = 0. ;
    return entry;
  }

  double b(int n,int m){
    double entry;
    if( 0 <= m && m <= n ){
      entry = sqrt( (n-m-1.)*(n-m) / ((2.*n-1.)*(2.*n+1.)) );
    } else if( -n <= m && m < 0 ){
      entry = -sqrt( (n-m-1.)*(n-m) / ((2.*n-1.)*(2.*n+1.)) );
    } else {
      entry = 0;
    }
    
    return entry;
  }
  

  //TODO : should be sttaic
  // Goes from (l,m) to (l',m')
  void Z_transfer::Compute(SpMatrix& coeff2, int L, double r, complex k, bool reg){
    double theta = 0;
    double phi = 0;
    
    // Initialization
    //clock_t begin = clock();
    complex val;
    for( int l = 0; l <= (2*L+1); l++ ){       
      int idx = index->LocateIndex(l,0); 
      
      complex F;
      (reg) ? F = Amos_sf_bessel_jl(l, k*r) : F = Amos_sf_hankel_1(l, k*r);
      val = sqrt(4*PI) * pow(-1., (double) l) * F * gsl_sf_harmonic(l,0,theta,phi);
      coeff2.add(idx, 0, val);
    }
    
    for( int n = 0; n <= (2*L+1); n++ ){       
      int idx = index->LocateIndex(n,0); 
      
      complex F;
      (reg) ? F = Amos_sf_bessel_jl(n, k*r) : F = Amos_sf_hankel_1(n, k*r);
      val = sqrt(4*PI) * F * gsl_sf_harmonic(n,0,theta,phi);
      coeff2.add(0, idx, val);
    }
    //cout << "Z_transfer initialization timing : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    
    
    // Compute sectorial ( |m| == n ) coefficients
    //TODO : use symmetry
    //begin = clock();
    for( int m = 0; m <= L; m++ ){
      for( int l = 1; l <= (2*L-m); l++ ){
	
	val = b(l,-m-1)/b(m+1,-m-1) * coeff2(index->LocateIndex(l-1,m), index->LocateIndex(m,m)) - 
	  b(l+1,m)/b(m+1,-m-1) * coeff2(index->LocateIndex(l+1,m), index->LocateIndex(m,m));
	coeff2.add(index->LocateIndex(l,m+1), index->LocateIndex(m+1,m+1), val);
	
	//TODO : is this if statement needed
	if( (l-1) >= m ){
	  val = b(l,-m-1)/b(m+1,-m-1) * coeff2(index->LocateIndex(l-1,-m), index->LocateIndex(m,-m)) - 
	    b(l+1,m)/b(m+1,-m-1) * coeff2(index->LocateIndex(l+1,-m), index->LocateIndex(m,-m));
	  coeff2.add(index->LocateIndex(l,-m-1), index->LocateIndex(m+1,-m-1), val);
	}
	
      }
    } 
    
    
    // Compute sectorial ( |s| == l ) coefficients
    for( int s = 0; s <= L; s++ ){
      for( int n = 1; n <= (2*L-s); n++ ){
	
	// TODO : is this if statement needed?
	if( (n-1) >= s ){
	  val = b(n,-s-1)/b(s+1,-s-1) * coeff2(index->LocateIndex(s,-s), index->LocateIndex(n-1,-s)) - 
	    b(n+1,s)/b(s+1,-s-1) * coeff2(index->LocateIndex(s,-s), index->LocateIndex(n+1,-s));
	  coeff2.add(index->LocateIndex(s+1,-s-1), index->LocateIndex(n,-s-1), val);
	}
	
	
	val = b(n,-s-1)/b(s+1,-s-1) * coeff2(index->LocateIndex(s,s), index->LocateIndex(n-1,s)) - 
	  b(n+1,s)/b(s+1,-s-1) * coeff2(index->LocateIndex(s,s), index->LocateIndex(n+1,s));
	coeff2.add(index->LocateIndex(s+1,s+1), index->LocateIndex(n,s+1), val);
	
      }
    }
    //cout << "Z_transfer sectorial timing : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    
    
    
    // Computation of remaining coefficients
    //begin = clock();
    for( int s = 0; s <= L; s++ ){
      int m = s;
      
      // Move forward in n
      for( int n = std::abs(m); n <= L; n++ ){
	for( int l = std::abs(s)+n-std::abs(m); l <= (2*L-n); l++ ){
	  val =  -a(l,s)/a(n,m) * coeff2(index->LocateIndex(l+1,s), index->LocateIndex(n,m));
	  
	  
	  if( (n-1) >= std::abs(m) )
	    val += a(n-1,m)/a(n,m) * coeff2(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	  
	  
	  if( (l-1) >= std::abs(s) )
	    val +=  a(l-1,s)/a(n,m) * coeff2(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	  
	  
	  coeff2.add(index->LocateIndex(l,s), index->LocateIndex(n+1,m), val);
	  coeff2.add(index->LocateIndex(n+1,-s), index->LocateIndex(l,-m), pow(-1., n+l+1.)*val);  // Symmetry
	  
	}
      }
      
      
      //Move forward in l
      for( int l = std::abs(s); l <= (L-std::abs(m)+std::abs(s)); l++ ){
	for( int n = std::abs(m)+l-std::abs(s); n <= 2*L-std::abs(m)-l+std::abs(s); n++ ){
	  
	  val = -a(n,m)/a(l,s) * coeff2(index->LocateIndex(l,s), index->LocateIndex(n+1,m));
	  
	  
	  if( (n-1) >= std::abs(m) )
	    val += a(n-1,m)/a(l,s) * coeff2(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	  
	  
	  if( (l-1) >= std::abs(s) )
	    val +=  a(l-1,s)/a(l,s) * coeff2(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	  
	  coeff2.add(index->LocateIndex(l+1,s), index->LocateIndex(n,m), val);
	  coeff2.add(index->LocateIndex(n,-s), index->LocateIndex(l+1,-m), pow(-1., n+l+1.)*val); // Symmetry
	  
	}
      }
      
    }
    //cout << "Z_transfer remaining timing : " << (clock() - begin)*1./CLOCKS_PER_SEC << endl;
    

    return;
  }
  
  
  // Proceed to transfer in-place
  void Z_transfer::Apply( std::vector<complex>& vec, Type type ) {
    
    int idx1, idx2;
    std::vector<complex> V;
    int k = 0;
    if( type == FORWARD ){
      assert( (int) vec.size() == (int) index2->size() );
      V.resize( index1->size() );
      
      for( int m = -L; m <= L; m++ ){
	for( int l1 = std::abs(m); l1 <= L_1; l1++ ){
	  for( int l2 = std::abs(m); l2 <= L_2; l2++ ){
	    idx1 = index1->LocateIndex(l1, m);
	    idx2 = index2->LocateIndex(l2, m);
	    
	    V[idx1] += coeff[k] * vec[idx2];
	    k++;
	  }
	}
      }
    } else if( type == ADJOINT ){
      assert( (int) vec.size() == (int) index1->size() );
      V.resize( index2->size() );
      
      for( int m = -L; m <= L; m++ ){
	for( int l1 = std::abs(m); l1 <= L_1; l1++ ){
	  for( int l2 = std::abs(m); l2 <= L_2; l2++ ){
	    idx1 = index1->LocateIndex(l1, m);
	    idx2 = index2->LocateIndex(l2, m);
	    
	    V[idx2] += conj(coeff[k]) * vec[idx1];
	    k++;
	  }
	}
      }
    }
    
    vec = V;
  }

  Indexing* Z_transfer::index;
  bool Z_transfer::Init = false;




  //---- Point & Shoot ----//
  // (Accelerated Translation/Transfer from 2 to 1)
  
  // Constructor
  PointAndShoot::PointAndShoot(double r, complex k_out, double theta, double phi,
			       Indexing* index1_, Indexing* index2_, bool reg)
    : index1(index1_), index2(index2_)
  {
    
    (index1->size() > index2->size()) ? Z_phi = new Z_Rotation(index1, phi) : Z_phi = new Z_Rotation(index2, phi);
    (index1->size() > index2->size()) ? Z_theta = new Z_Rotation(index1, theta) : Z_theta = new Z_Rotation(index2, theta);
    
    S_F = new S_Rotation();
    S_B = new S_Rotation();
    
    Z_T = new Z_transfer(r, k_out, index1, index2, reg);  
  }
   

  // Destructor 
  PointAndShoot::~PointAndShoot()
  {
    delete Z_phi;
    delete Z_theta;
    
    delete S_F;
    delete S_B;
    
    delete Z_T;
  }
  
  
  // Note : pass by value since rotations act in place
  std::vector<complex> PointAndShoot::Apply( std::vector<complex>& v, Type type )
  {
    // TODO: HACK TO AVOID CHANGING SIZE OF V
    std::vector<complex> vec = v;

    ( type == FORWARD ) ? assert( (int) vec.size() == (int) index2->size() ) : assert( (int) vec.size() == (int) index1->size() );  
    int L;


    //--- Forward rotation ---// 
    ( type == FORWARD ) ? index = index2 : index = index1;
    ( type == FORWARD ) ? L = (*index2)(index2->size()-1,0) :  L = (*index1)(index1->size()-1,0);
    
    //Perform positive phi rotation
    Z_phi->Apply(vec, index, Z_Rotation::POS);

    //Apply S matrix
    S_F->Apply(vec, L, FORWARD);
    
    //Perform positive theta rotation
    Z_theta->Apply(vec, index, Z_Rotation::POS);
    
    //Apply Hermitian transpose of S matrix
    S_F->Apply(vec, L, ADJOINT);
    
    //---Transfer---//
    //Apply transfer along the z-axis
    Z_T->Apply(vec, type);
    
    ( type == FORWARD ) ? assert( (int) vec.size() == (int) index1->size() ) : assert( (int) vec.size() == (int) index2->size() );  
    
    //---Backward rotation---//
    ( type == FORWARD ) ? index = index1 : index = index2;
    ( type == FORWARD ) ? L = (*index1)(index1->size()-1,0) :  L = (*index2)(index2->size()-1,0);
    
    //Apply S matrix
    S_B->Apply(vec, L, FORWARD);
    
    //Perform negative theta rotation
    Z_theta->Apply(vec, index, Z_Rotation::NEG);
    
    //Apply Hermitian transpose of S matrix
    S_B->Apply(vec, L, ADJOINT);
    
    //Perform negative phi rotation
    Z_phi->Apply(vec, index, Z_Rotation::NEG);
    
    return vec;
  }
  
  



  
  
  // ***Full translation matrix computation
  // Computes full transfer matrix from 2 to 1 along the std::vector (r,theta,phi)
  
  Gumerov_TransCoeff::Gumerov_TransCoeff(double r, complex k_out, double theta, double phi, int L_, bool reg) : L(L_)
  {
    TCMatrix = Gumerov_TransCoeff::Compute(L, r, theta, phi, k_out, reg);
    
    // Output to file
    if( (theta - PI/2) < 0 ) {
      ofstream file_r("TransferMatrix_real_r=1.csv", ios::out);
      ofstream file_i("TransferMatrix_imag_r=1.csv", ios::out);
      //cout << "[" ;
      for( int i = 0; i < (L+1)*(L+1); i++ ){
	for( int j = 0; j < (L+1)*(L+1); j++ ) {
	  //if( j != 0 ) { cout << ","; }
	  if( j != 0 ) { file_r << ","; }
	  if( j != 0 ) { file_i << ","; }
	  // cout << real( TCMatrix(i,j) ) << "+" <<
	  //	imag( TCMatrix(i,j) )  << "*1i";
	  file_r << std::setprecision(15) << real( TCMatrix(i,j) ) ;
	  file_i << std::setprecision(15) << imag( TCMatrix(i,j) ) ;
	}
	//cout << ";" << endl;
	file_r << endl;
	file_i << endl;
      }
      //cout << "]" << endl << endl;
      
      file_r.close();
      file_i.close();
    }
  }
  
  /*std::vector<complex>  Gumerov_TransCoeff::Apply( std::vector<complex> vec, Type type = FORWARD )
  {
    std::vector<complex> vec_out(vec.size());
    
    if( type == FORWARD ){
      for( int i = 0; i < (L+1)*(L+1); i++ )
	for( int j = 0; j < (L+1)*(L+1); j++ )
	  vec_out[i] += TCMatrix(i,j)*vec[j];
      
    } else if( type == ADJOINT ){
      for( int i = 0; i < (L+1)*(L+1); i++ ){
	for( int j = 0; j < (L+1)*(L+1); j++ ){
	  vec_out[i] += conj(TCMatrix(j,i))*vec[j];
	}
      }
      
    }
    
    return vec_out;
  }*/
  
  
  
  
  
  
  // Goes from (l,m) to (l',m')
  arma::Mat<complex>  Gumerov_TransCoeff::Compute(int L, double r, double theta, double phi, complex k, bool reg){
    arma::Mat<complex> M((2*L+5)*(2*L+5), (2*L+5)*(2*L+5));
    Indexing *index;
    
    // Could have a static index and update only when needed. Update only requires adding elements
    index = new Indexing(2*L+5);
    
    // Initialization
    complex val;
    for( int l = 0; l <= (2*L+1); l++ ){   
      for( int s = -l; s <= l; s++ ){
	complex F;
	(reg) ? F = Amos_sf_bessel_jl(l, k*r) : F = Amos_sf_hankel_1(l, k*r);
	val = sqrt(4*PI) * pow(-1., (double) l) * F * gsl_sf_harmonic(l,-s,theta,phi);
	
	( -s > 0 ) ? M(index->LocateIndex(l,s), 0) = val*pow(-1., (double) s) : M(index->LocateIndex(l,s), 0) = val;
      }
    }
    
    for( int n = 0; n <= (2*L+1); n++ ){       
      for( int m = -n; m <= n; m++ ){
	complex F;
	(reg) ? F = Amos_sf_bessel_jl(n, k*r) : F = Amos_sf_hankel_1(n, k*r);
	val = sqrt(4*PI) * F * gsl_sf_harmonic(n,m,theta,phi);
	
	( m > 0 ) ? M(0, index->LocateIndex(n,m)) = val*pow(-1., (double) m) : M(0, index->LocateIndex(n,m)) = val;
      }
    }
    
    
    // Compute sectorial ( |m| == n ) coefficients
    for( int m = 0; m <= L; m++ ){
      for( int l = 0; l <= (2*L-m); l++ ){
	for( int s = -l; s <= l; s++ ){
	  
	  
	  if( !( l==0 && s==0 ) ){
	    if( l-1 >= std::abs(s-1) ) 
	      val = b(l,-s)/b(m+1,-m-1) * M(index->LocateIndex(l-1,s-1), index->LocateIndex(m,m))
		- b(l+1,s-1)/b(m+1,-m-1) * M(index->LocateIndex(l+1,s-1), index->LocateIndex(m,m));
	    else
	      val = - b(l+1,s-1)/b(m+1,-m-1) * M(index->LocateIndex(l+1,s-1), index->LocateIndex(m,m));
	    
	    M(index->LocateIndex(l,s), index->LocateIndex(m+1,m+1)) = val;
	    
	    
	    if( l-1 >= std::abs(s+1) ) 
	      val = b(l,s)/b(m+1,-m-1) * M(index->LocateIndex(l-1,s+1), index->LocateIndex(m,-m))
		-b(l+1,-s-1)/b(m+1,-m-1) * M(index->LocateIndex(l+1,s+1), index->LocateIndex(m,-m));
	    else
	      val = -b(l+1,-s-1)/b(m+1,-m-1) * M(index->LocateIndex(l+1,s+1), index->LocateIndex(m,-m));
	    
	    M(index->LocateIndex(l,s), index->LocateIndex(m+1,-m-1)) = val;
	    
	  }
	  
	}
      }
    }
    
    
    // Compute sectorial ( |s| == l ) coefficients
    for( int s = 0; s <= L; s++ ){
      for( int n = 0; n <= (2*L-s); n++ ){
	for( int m = -n; m <= n; m++ ){
	  
	  if( !( n==0 && m==0 ) ){
	    
	    if( n-1 >= std::abs(m+1) ) 
	      val = b(n,m)/b(s+1,-s-1) * M(index->LocateIndex(s,-s), index->LocateIndex(n-1,m+1))
		- b(n+1,-m-1)/b(s+1,-s-1) * M(index->LocateIndex(s,-s), index->LocateIndex(n+1,m+1));
	    else
	      val = - b(n+1,-m-1)/b(s+1,-s-1) * M(index->LocateIndex(s,-s), index->LocateIndex(n+1,m+1));
	    
	    M(index->LocateIndex(s+1,-s-1), index->LocateIndex(n,m)) = val;
	    
	    
	    if( n-1 >= std::abs(m-1) ) 
	      val = b(n,-m)/b(s+1,-s-1) * M(index->LocateIndex(s,s), index->LocateIndex(n-1,m-1))
		-b(n+1,m-1)/b(s+1,-s-1) * M(index->LocateIndex(s,s), index->LocateIndex(n+1,m-1));
	    else
	      val = -b(n+1,m-1)/b(s+1,-s-1) * M(index->LocateIndex(s,s), index->LocateIndex(n+1,m-1));
	    
	    M(index->LocateIndex(s+1,s+1), index->LocateIndex(n,m)) = val;
	    
	  }
	  
	} 
      }
    }
    
    // Computation of remaining coefficients
    for( int s = -L; s <= L; s++ ){
      for( int m = -L; m <= L; m++ ){
	
	if( std::abs(s) <= std::abs(m) ){
	  
	  // Move forward in n
	  for( int n = std::abs(m); n <= L; n++ ){
	    for( int l = std::abs(s)+n-std::abs(m); l <= (2*L-n); l++ ){
	      
	      val =  -a(l,s)/a(n,m) * M(index->LocateIndex(l+1,s), index->LocateIndex(n,m));
	      
	      if( (n-1) >= std::abs(m) )
		val += a(n-1,m)/a(n,m) * M(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	      
	      if( (l-1) >= std::abs(s) )
		val +=  a(l-1,s)/a(n,m) * M(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	      
	      
	      M(index->LocateIndex(l,s), index->LocateIndex(n+1,m)) =  val;
	      
	      
	    }
	  }
	  
	  
	  //Move forward in l
	  for( int l = std::abs(s); l <= (L-std::abs(m)+std::abs(s)); l++ ){
	    for( int n = std::abs(m)+l-std::abs(s); n <= 2*L-std::abs(m)-l+std::abs(s); n++ ){
	      
	      val = -a(n,m)/a(l,s) * M(index->LocateIndex(l,s), index->LocateIndex(n+1,m));
	      
	      if( (n-1) >= std::abs(m) )
		val += a(n-1,m)/a(l,s) * M(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	      
	      
	      if( (l-1) >= std::abs(s) )
		val +=  a(l-1,s)/a(l,s) * M(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	      
	      M(index->LocateIndex(l+1,s), index->LocateIndex(n,m)) = val;
	    }
	  }
	  
	} else if( std::abs(s) > std::abs(m) ){
	  
	  
	  //Move forward in l
	  for( int l = std::abs(s); l <= L; l++ ){
	    for( int n = std::abs(m)+l-std::abs(s); n <= 2*L-1; n++ ){
	      
	      val = -a(n,m)/a(l,s) * M(index->LocateIndex(l,s), index->LocateIndex(n+1,m));
	      
	      if( (n-1) >= std::abs(m) )
		val += a(n-1,m)/a(l,s) * M(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	      
	      if( (l-1) >= std::abs(s) )
		val +=  a(l-1,s)/a(l,s) * M(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	      
	      M(index->LocateIndex(l+1,s), index->LocateIndex(n,m)) = val;
	    }
	  }
	  
	  // Move forward in n
	  for( int n = std::abs(m); n <= L-std::abs(s)+std::abs(m) ; n++ ){
	    for( int l = std::abs(s)+n-std::abs(m); l <= (2*L-std::abs(s)-n+std::abs(m)); l++ ){
	      val =  -a(l,s)/a(n,m) * M(index->LocateIndex(l+1,s), index->LocateIndex(n,m));
	      
	      if( (n-1) >= std::abs(m) )
		val += a(n-1,m)/a(n,m) * M(index->LocateIndex(l,s), index->LocateIndex(n-1,m));
	      
	      if( (l-1) >= std::abs(s) )
		val +=  a(l-1,s)/a(n,m) * M(index->LocateIndex(l-1,s), index->LocateIndex(n,m));
	      
	      
	      M(index->LocateIndex(l,s), index->LocateIndex(n+1,m)) =  val;
	    }
	  }
	  
	  
	}
	
	
      }	
    }
    
    int N = (L+1)*(L+1);
    M.resize( N, N); // Only keep desired entries
    
    // Normalization required due to discrepancy between current convention for
    // spherical harmonics and Gumerov's convention
    for( int i = 0; i < N; i++ ){
      for( int j = 0; j < N; j++ ){
	double c = 1.;
	if( (*index)(i,1) >= 0 ) { c *= pow(-1., (double) (*index)(i,1)); }
	if( (*index)(j,1) >= 0 ) { c *= pow(-1., (double) (*index)(j,1)); }
	
	M(i,j) *= c;
	//cout << "," << std::setprecision(15) <<  M(i,j);
      }
      //cout << ";" << endl;
    }
    //cout << endl << endl;
    
    delete index;
    return M;
  }
  
  
   


 


  /*
// r: radius of scatterers
// R: 
  int L_MAX_COMPUTE(complex k_out, complex k, double d){

    //srand48(time(NULL));

    int N = 30;
    double R = d;
    double theta = PI/4.;
    double phi = PI/4.;
   
    Gumerov_TransCoeff TransMat(R, k_out, theta, phi, N, false);
    arma::Mat<complex> Mat = TransMat.TCMatrix;

    // Apply T-matrix
    for( int l_p = 0; l_p <= N; l_p++ ){
      for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	int row_idx = l_p*l_p + (l_p+m_p);
	
	for( int l = 0; l <= N; l++ ){
	  for( int m = -l; m <= l; m++ ){
	    complex mie = Mie(k, k_out, RADIUS, RHO, l_p, m_p) * Mie(k, k_out, RADIUS, RHO, l, m) / Amos_sf_bessel_jl(l_p, k_out*d/4.);
	    int col_idx = l*l + (l+m);
	    
	    Mat(row_idx, col_idx) *= mie;
	  }
	}
      }
    }
    double norm = arma::norm(Mat, 2);
    cout << "norm :: " << norm << endl;
    int L = 2;
    //int rank = arma::rank(Mat, EPS); 

      int L = 1;//stc::ceil(sqrt(rank));
    // cout << "rank : " << rank << endl;
    while( L <= N ){
      Gumerov_TransCoeff TransMat2(R, K_OUT, theta, phi, L, false);
      arma::Mat<complex> Mat2 = TransMat2.TCMatrix;
      
      // Apply T-matrix
      for( int l_p = 0; l_p <= L; l_p++ ){
	for( int m_p = -l_p; m_p <= l_p; m_p++ ){
	  int row_idx = l_p*l_p + (l_p+m_p);
	  complex mie = Mie(K, K_OUT, RADIUS, RHO, l_p, m_p);
	  
	  for( int l = 0; l <= L; l++ ){
	    for( int m = -l; m <= l; m++ ){
	      int col_idx = l*l + (l+m);
	      
	      Mat2(row_idx, col_idx) *= mie;
	    }
	  }
	}
      }

      int P = arma::rank(Mat2, EPS);
      cout << " P : " << P << endl;
      if( P == (L+1)*(L+1) )
	L++;
      else
	break;
    }

    //cout << L << " = ? " <<  std::ceil(sqrt( arma::rank(Mat2, EPS) )) << endl;
    

    cout << "NEW L : " <<  L << endl;
    return L;
  }


}

  */

  











}

#endif
