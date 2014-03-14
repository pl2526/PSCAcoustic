#ifndef SOLVER_CPP
#define SOLVER_CPP

/*
 *  Solver.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 3/6/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Iterative linear solver implementation 
 *  through PETSc library.
 */


#include "Solver.h"

// TODO: Should not have different instances of Index


// Compute L^2 norm
double norm(std::vector< std::vector<complex > >& vec)
{

  double val = 0;
  for( int i = 0; i < (int) vec.size(); i++ )
    for( int j = 0; j < (int) vec[i].size(); j++ )
      val += std::abs(vec[i][j])*std::abs(vec[i][j]);
  
  return sqrt(val);
}


Solver::Solver(){ Index.Compute(LMAX); }


void Solver::RHSInit( IncWave* IW, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& RHS){
  std::vector<IncWave*> vec_IW(1); vec_IW[0] = IW; 
  RHSInit(  vec_IW, ScL, RHS);
}


//The right hand side of the system corresponds to the incostd::ming wave multipole
//coefficients for an expansion about the center of the current scatterer
void Solver::RHSInit( std::vector<IncWave*> IW, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& RHS)
{

  // Initialize
  for( int n = 0 ; n < NScat ; n++ )
    for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ )
      RHS[n][i] = 0;
 

  int N_src = IW.size();
  // Scratch space
  std::vector<complex> rhs(RHS[0].size());
  // for each source
  for( int k = 0; k < N_src; k++ ){
    
    // For each scatterer
    for( int n = 0 ; n < NScat ; n++ )
      {
	//Transfer expansion about center of scatterer
	Pvec p = ScL[n].getLocation();
	IW[k]->Translate(p, rhs);

	//Apply T-matrix
	ScL[n].TM_Apply(rhs);

	for( int i = 0; i < (LMAX+1)*(LMAX+1); i++ )
	  RHS[n][i] += rhs[i];
	
      }
  }
  
  return;
}


// Solve linear solver, variant of GMRES
void Solver::Solve(PSC_Env& PSC_env, std::vector<Scatterer>& ScL, std::vector< std::vector<complex> >& RHS, 
		   std::vector< std::vector<complex> >& u,  double& res, double& rel_res, double& cond )
{
    //u = RHS;
    //return;
    

   // Right-hand side
    Matrix<complex, Dynamic, 1> b; b.resize(NScat*Index.size(),1);
    std::vector< std::vector<complex> > v(NScat, std::vector<complex>(Index.size()));
    Matrix<complex, Dynamic, 1> x;
    Matrix<complex, Dynamic, 1> y;
    Matrix<complex, Dynamic, Dynamic> H; H.resize(2,1); // Hessenberg matrix    

    double norm_RHS = norm(RHS);
    for( int n = 0; n < NScat; n++ ){
      for( int i = 0; i < Index.size(); i++ ){
	b(n*Index.size()+i) = RHS[n][i];
	u[n][i] = RHS[n][i];///norm_RHS;
      }
    }



  if(0){
    // *** GMRES ***
    Matrix<complex, Dynamic, Dynamic> R; R.resize(NScat*Index.size(),1);
    Matrix<complex, Dynamic, Dynamic> Q; Q.resize(NScat*Index.size(),1);
    //FullPivHouseholderQR<complex, Dynamic, Dynamic> QR;
    Matrix<complex, Dynamic, 1> z; 

    // Initialization
    // Apply scattering operator
    cout << "    ***Applying scattering operator (TA)..." << endl;
    ScatteringOperator(PSC_env, ScL, u, v);
    cout << "    ***Applying scattering operator: done" << endl;
    double rho_0 = 0;
    for( int n = 0; n < NScat; n++ ){
      for( int i = 0; i < Index.size(); i++ ){
	Q(n*Index.size() + i) = u[n][i] - v[n][i];
	rho_0 += std::abs(Q(n*Index.size() + i)) * std::abs(Q(n*Index.size() + i)) ;
      }
    }
    rho_0 = sqrt(rho_0);
    for( int n = 0; n < NScat; n++ )
      for( int i = 0; i < Index.size(); i++ )
	Q(n*Index.size() + i) /= rho_0;


    double res = 1e10;
    int k = 0;
    while( k < MAXITS ){
      if( k > 0 ){
	for( int n = 0; n < NScat; n++ )
	  for( int i = 0; i < Index.size(); i++ )
	    Q(n*Index.size() + i,k) = R(n*Index.size() + i,k-1) / H(k,k-1);
      }

      for( int n = 0; n < NScat; n++ ){
	for( int i = 0; i < Index.size(); i++ ){
	  u[n][i] = Q(n*Index.size() + i,k);
	  v[n][i] = 0.;
	}
      }
      cout << "    ***Applying scattering operator (TA)..." << endl;
      ScatteringOperator(PSC_env, ScL, u, v);
      cout << "    ***Applying scattering operator: done" << endl;
      for( int n = 0; n < NScat; n++ )
	for( int i = 0; i < Index.size(); i++ )
	  R(n*Index.size() + i,k) = u[n][i] - v[n][i];
      
      // Orthogonalize
      for( int j = 0; j <= k; j++ ){
	H(j,k) = 0;
	  for( int n = 0; n < NScat; n++ )
	    for( int i = 0; i < Index.size(); i++ )
	      H(j,k) += std::conj(Q(n*Index.size() + i,j)) * R(n*Index.size() + i,k);

	  for( int n = 0; n < NScat; n++ )
	    for( int i = 0; i < Index.size(); i++ )
	      R(n*Index.size() + i,k) -= H(j,k) * Q(n*Index.size() + i,j);
      }
      H(k+1,k) = 0;
      complex val = 0;
      for( int n = 0; n < NScat; n++ ){
	for( int i = 0; i < Index.size(); i++ ){
	  H(k+1,k) += std::abs(R(n*Index.size() + i,k)) * std::abs(R(n*Index.size() + i,k));
	  if( k > 1 )
	    val += std::conj(Q(n*Index.size() + i,0)) * Q(n*Index.size() + i,k);
	}
      }
      H(k+1,k) = sqrt(H(k+1,k));
      
      cout << "val : " << val << endl;

      /* cout << "before" << endl;
      	cout << H << endl;
	cout << endl << endl << H.rows() << " : " << H.cols() << endl;*/
      
      // Compute L-S residual
      if( k%K_LS == 0 ){

	// R.h.s.
	z.resize(H.rows(),1);
	for( int i = 0; i < H.rows(); i++ )
	  z(i) = 0.;
	z(0) = rho_0;

	// Solve L-S
	y = H.fullPivHouseholderQr().solve(z);

	Matrix<complex,Dynamic,1> nu = H*y - z;
	Matrix<complex,1,1> mu = nu.adjoint() * nu;
	res = std::abs(sqrt(mu(0)));
	rel_res = res / norm_RHS;
	
	cout << "L-S residual : " << std::setprecision(15) << res << endl;
	cout << "L-S relative residual : " << std::setprecision(15) << rel_res << endl << endl;
	
	if( res < ATOL || rel_res < RTOL){
	  cout << "Stopping criterion satisfied" << endl;
	  
	  x = b + Q*y;
	  for( int n = 0; n < NScat; n++ )
	    for( int i = 0; i < (int) Index.size(); i++ )
	      u[n][i] += norm_RHS*x(n*Index.size() + i);

	  //H.topLeftCorner(H.cols(),H.cols());

	  
	  break;
	}
      }

      
      
      Q.conservativeResize(NScat*Index.size(), Q.cols()+1);
      R.conservativeResize(NScat*Index.size(), R.cols()+1);
      H.conservativeResize(H.rows()+1, H.cols()+1);

      k++;
    }


  } else {

 
    
    

    
    Matrix<complex, Dynamic, Dynamic> A; A.resize(NScat*Index.size(),1);
    Matrix<complex, Dynamic, Dynamic> M; M.resize(NScat*Index.size(),1);

    
    // Solution
    int k = 0;
    double comp_time = clock();
    while( k < MAXITS )
      {
	cout << "Iteration : " << k << endl;
	
	
	// Hessenberg matrix
	// Put zero values to zero in Hessenberg matrix (for computing condition number estimate)
	//H.conservativeResize((k+2), k+1);
	//for( int j = 0; j < (k-1); j++ )
	// H(k,j) = 0;
	
	// Orthogonalize and apply operator
	for( int j = 0; j < k; j++ )
	  {
	    complex dot = 0.;
	    for( int n = 0; n < NScat; n++ ){
	      for( int i = 0; i < Index.size(); i++ ){
		dot += conj( M(n*Index.size() + i, j) ) * u[n][i];
	      }
	    }
	    // Populate Hessenberg matrix
	    //H(j,k) = dot;
	    
	    for( int n = 0; n < NScat; n++ ){
	      for( int i = 0; i < Index.size(); i++ ){
		u[n][i] -= dot*M(n*Index.size() + i, j);
	      }
	    }

	  }
	
	// Normalize
	double Norm = norm(u);
	//H(k+1,k) = Norm;
	for( int n = 0; n < NScat; n++ ){
	  for( int i = 0; i < Index.size(); i++ ){
	    u[n][i] /= Norm;
	    M(n*Index.size()+i, k) = u[n][i];
	  }
	}


	for( int j = 0; j < k; j++ ){
	  complex dot = 0;
	  for( int n = 0; n < NScat; n++ )
	    for( int i = 0; i < Index.size(); i++ )
	      dot += conj( M(n*Index.size() + i, j) ) * M(n*Index.size() + i, k);
	}
	      

	
	// Apply scattering operator
	cout << "    ***Applying scattering operator (TA)..." << endl;
	ScatteringOperator(PSC_env, ScL, u, v);
	cout << "    ***Applying scattering operator: done" << endl;
	//cout << "norm : " << norm(v) << endl;



	// TODO: could be faster
	// Store values and build linear system
	for( int n = 0; n < NScat; n++ )
	  for( int i = 0; i < Index.size(); i++ )
	    A(n*Index.size()+i, k) = u[n][i] - v[n][i];


	
	// Update
	u = v;


	// Compute L-S residual
	if( k%K_LS == 0 ){

	  // Solve L-S
	  x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

	  y = A*x-b;
	  Matrix<complex,1,1> nu = y.adjoint() * y;
	  res = abs(sqrt(nu(0,0)));
	  rel_res = res / norm(RHS);

	  cout << "L-S residual : " << std::setprecision(15) << res << endl;
	  cout << "L-S relative residual : " << std::setprecision(15) << rel_res << endl << endl;

	  if( res < ATOL || rel_res < RTOL){
	    cout << "Stopping criterion satisfied" << endl;

	    // TODO: currently inefficient ?
	    for( int n = 0; n < (int) NScat; n++ )
	      for( int i = 0; i < (int) Index.size(); i++ )
		u[n][i] = 0.;
	        
	    for( int k = 0; k < M.cols(); k++)
	      for( int n = 0; n < NScat; n++ )
		for( int i = 0; i < (int) Index.size(); i++ )
		  u[n][i] += x(k) * M(n*Index.size()+i, k);
	    
	    break;
	  }
	}

	
	// Increase size if necessary
	A.conservativeResize(NScat*Index.size(), A.cols()+1);
	M.conservativeResize(NScat*Index.size(), M.cols()+1);
	
	k++;    // Increment iterate
      }
    comp_time = (clock() - comp_time)*1./CLOCKS_PER_SEC;
    
    
    
    // Check accuracy
    double acc = 0;
    ScatteringOperator(PSC_env, ScL, u, v);
    for( int n = 0; n < (int) NScat; n++ )
      for( int i = 0; i < (int) Index.size(); i++ )
	acc += pow(std::abs( u[n][i] - v[n][i] - RHS[n][i] ), 2.);
    cout << "Accuracy : " << sqrt(acc) << endl;
    cout << "Relative accuracy : " << sqrt(acc)/norm(u) << endl;
    
    
    // Compute estimate of condition number
    JacobiSVD< Matrix<complex, Dynamic, Dynamic> > svd(A, ComputeThinU | ComputeThinV);
    //JacobiSVD< Matrix<complex, Dynamic, Dynamic> > svd(H.rightCols(H.cols()-1), ComputeThinU | ComputeThinV);
    //Eigen::SelfAdjointEigenSolver< Matrix<complex, Dynamic, Dynamic> > svd(H);
    
    cout << "FOO" << endl;
    
    double sigma_max = 0;
    double sigma_min = 1e10;
    for( int i = 0; i < (A.cols()-1); i++ ){
      sigma_max = std::max(sigma_max, std::abs(1.-(svd.singularValues())(i)) );
      sigma_min = std::min(sigma_min, std::abs(1.-(svd.singularValues())(i)) );
    }
    
    cond = sigma_max / sigma_min;
    cout << "Approximation of condition number : " << cond << endl;
    
  }
}

// Operator TA
void Solver::ScatteringOperator( PSC_Env& PSC_env, std::vector<Scatterer>& ScL, std::vector< std::vector<complex > >& u,
				 std::vector< std::vector<complex > >& v){

  for( int n = 0; n < NScat; n++ )
    for( int i = 0; i < Index.size(); i++ )
      v[n][i] = 0;

  // Compute interactions among scatterers for each type of waves
  PSC_env.getTransfer()->execute(u, v);

  // Apply T-matrix
  for( int n = 0; n < NScat ; n++ )
    ScL[n].TM_Apply(v[n]);

  return;
}







#endif