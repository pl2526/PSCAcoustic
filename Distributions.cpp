/*
 *  Districutions.h
 *  PSCAcoustic
 *
 *  Various routines to generate arrangements of scatterers.
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


#ifndef DISTRIBUTIONS_CPP
#define DISTRIBUTIONS_CPP

#include "General.h"
#include "Distributions.h"


//Return std::vector of scatterers with location and properties assigned below
std::vector<Scatterer> ScatInit( int N, complex k, complex k_out ){
  
  std::vector<Scatterer> ScL;
  
  if( 0 ){

    // Compute center locations
    Pvec CENTER(0.,0.,0. ,Pvec::Cartesian);                                          // Center of cluster
    std::vector<Pvec> center = RandSpherical(0.3, RADIUS, D_MIN, CENTER, NScat);    // Random distribution
    
    assert( center.size() == NScat );

    // Construct scatterers
    for( int n = 0; n < NScat; n++ ){
      Scatterer scatterer(RADIUS, k, k_out, RHO, center[n]);
      //cout << center[n].x << " : " << center[n].y << " : " << center[n].z << " : " << endl;
      ScL.push_back(scatterer);
    }


  }else {
    
    //Initialize location of scatterers
    
    Pvec vec0(0., 0.0, 0.0, Pvec::Cartesian); 
    Scatterer scatterer0(RADIUS, k, k_out, RHO, vec0);
    ScL.push_back(scatterer0);
    
    /*Pvec vec1(0., -0.45, 0., Pvec::Cartesian); 
    Scatterer scatterer1(RADIUS, k, k_out, RHO, vec1);
    ScL.push_back(scatterer1);*/

    /*Pvec vec2(-3., -2.5, 0, Pvec::Cartesian); 
    Scatterer scatterer2(RADIUS, K, K_OUT, RHO, vec2);
    ScL.push_back(scatterer2);*/
  }
  
  cout << "NScat : " << ScL.size() << endl;
  assert( NScat == ScL.size() );
  return ScL;
}








//( Used only by ScatInit() )
// Returns a std::vector of location of sphere randomly located around the origin within
// a spherical enclosure
std::vector<Pvec> RandSpherical(double R, double r, double d, Pvec c, double N){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  
  //srand48(0);
  srand48(time(0));
  double dist = 0;
  double x,y,z;
  int maxN = 5000*N; //Maximum number of iterations

  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {
      x = 2*R*drand48() - R;
      y = 2*R*drand48() - R;
      z = 2*R*drand48() - R;
      dist = sqrt( x*x + y*y + z*z );
      //cout << "dist : " << dist << endl;
      
      if( (dist+r) <= R ){
	//Check sphere is far enough from others
	bool too_close = false;
	int i = 0;
	while( i < (int) array.size() && !too_close ){

	  dist = sqrt( (x-array[i].x)*(x-array[i].x) + 
		       (y-array[i].y)*(y-array[i].y) + 
		       (z-array[i].z)*(z-array[i].z) );

	  //cout << "(" << array[i].x << "," << array[i].y << "," << array[i].z << ")" << endl;
	  //cout << "dist2 : " << dist << endl;

	  if( dist < (2*r+d) )
	    too_close = true;

	  i++;
	}

	if( !too_close ){
	  //cout << "(" << x << "," << y << "," << z << ")" << endl;
	  Pvec center(x, y, z, Pvec::Cartesian);
	  array.push_back(c + center);
	}	

      }

      n++;
  }
  
  return array;
  
  }


std::vector<Pvec> TruncatedSphere(double R, double t, double r, double d, Pvec c, double N, Pvec src_loc, bool rand_flag){
  // R: radius of sphere
  // t: distance between center and cut (along y-axis) 
  // r: radius of the spheres
  // d: minimum distance allowable between any sphere

  if( rand_flag )
    srand48(time(0));
  else
    srand48(0);

  
  std::vector<Pvec> array; // Array of center locations
  double dist = 0;
  double dist_src = 0;
  double x,y,z;
  int maxN = 5000*N; //Maximum number of iterations

  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {
      x = 2*R*drand48() - R;
      y = 2*R*drand48() - R;
      z = 2*R*drand48() - R;
      dist = sqrt( x*x + y*y + z*z );
      

      
      if( (dist + r) <= R && y < t){
	
	//Check sphere is far enough from others and source
	bool too_close = false;
	int i = 0;
	while( i < (int) array.size() && !too_close ){

	  dist = sqrt( (x-array[i].x)*(x-array[i].x) + 
		       (y-array[i].y)*(y-array[i].y) + 
		       (z-array[i].z)*(z-array[i].z) );

	  // Check is sphere is far enough from source 
	  dist_src = sqrt( (x-src_loc.x)*(x-src_loc.x) + 
		       (y-src_loc.y)*(y-src_loc.y) + 
		       (z-src_loc.z)*(z-src_loc.z) );



	  //cout << "(" << array[i].x << "," << array[i].y << "," << array[i].z << ")" << endl;
	  //cout << "dist2 : " << dist << endl;

	  if( dist < (2*r+d) || dist_src < (4*r+2.*d))
	    too_close = true;

	  i++;
	}

	if( !too_close ){
	  //cout << "(" << x << "," << y << "," << z << ")" << endl;
	  Pvec center(x, y, z, Pvec::Cartesian);
	  array.push_back(c + center);
	}	

      }

      n++;
  }
  
  return array;
  
  }


std::vector<Pvec> ErgodicCavity(double R, double t, double thick, Pvec center, double r, double d, double N, Pvec src_loc, bool rand_flag){
  
  //if( rand_flag )
    //srand48(time(0));
  //else
    srand48(0);

  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  double dist_src = 0;
  int maxN = 100*N; //Maximum number of iterations
  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {

      std::vector<double> Z(3);
      Z[0] = (2*R*drand48() - R) + center.x;
      Z[1] = (2*R*drand48() - R) + center.y;

      // Spheroid
      //Z[2] = (2*R*drand48() - R) + center.z;

      // Ergodic cylinder
      Z[2] = (2*thick*drand48() - thick) + center.z;

      //Check sphere is far enough from others
      bool too_close = false;
      double dist;
      int i = 0;
      // Creates edge
      double rad = sqrt( (Z[0]-center.x)*(Z[0]-center.x) +
			 (Z[1]-center.y)*(Z[1]-center.y) );
      if( Z[1] < t && rad <= R ){
	while( i < (int) array.size() && !too_close ){
	  
	  dist = sqrt( (Z[0]-array[i].x)*(Z[0]-array[i].x) + 
		       (Z[1]-array[i].y)*(Z[1]-array[i].y) + 
		       (Z[2]-array[i].z)*(Z[2]-array[i].z) );

	  dist_src = sqrt( (Z[0]-src_loc.x)*(Z[0]-src_loc.x) + 
			   (Z[1]-src_loc.y)*(Z[1]-src_loc.y) + 
			   (Z[2]-src_loc.z)*(Z[2]-src_loc.z) );

	  
	  if( dist < (2*r+d) || dist_src < (2.*r+2.*d) )
	    too_close = true;
	  
	  i++;
	}
	
	if( !too_close ){
	  Pvec loc(Z[0], Z[1], Z[2], Pvec::Cartesian);
	  array.push_back(loc);
	}	
      }

      n++;
    }
  
  cout << "size : " << array.size() << endl;
  return array;  
}



std::vector<Pvec> RandSphericalXY(double R, double r, double d, Pvec c, double N){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  
  //srand48(0);
  srand48(time(0));
  double dist = 0;
  double x,y,z;
  int maxN = 100*N; //Maximum number of iterations

  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {
      x = 2*R*drand48() - R;
      y = 2*R*drand48() - R;
      z = 0.;
      dist = sqrt( x*x + y*y + z*z );
      //cout << "dist : " << dist << endl;
      
      if( (dist+r) <= R ){
	//Check sphere is far enough from others
	bool too_close = false;
	int i = 0;
	while( i < (int) array.size() && !too_close ){

	  dist = sqrt( (x-array[i].x)*(x-array[i].x) + 
		       (y-array[i].y)*(y-array[i].y) + 
		       (z-array[i].z)*(z-array[i].z) );

	  //cout << "(" << array[i].x << "," << array[i].y << "," << array[i].z << ")" << endl;
	  //cout << "dist2 : " << dist << endl;

	  if( dist < (2*r+d) )
	    too_close = true;

	  i++;
	}

	if( !too_close ){
	  //cout << "(" << x << "," << y << "," << z << ")" << endl;
	  Pvec center(x, y, z, Pvec::Cartesian);
	  array.push_back(c + center);
	}	

      }

      n++;
  }
  
  return array;
  
}


std::vector<Pvec> RandRectangularXY(double x_dim, double y_dim, double r, double d, Pvec c, double N){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  
  //srand48(0);
  srand48(time(0));
  double dist = 0;
  double x,y,z;
  int maxN = 100*N; //Maximum number of iterations

  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {
      x = 2*x_dim*drand48() - x_dim;
      y = 2*y_dim*drand48() - y_dim;
      z = 0.;
      dist = sqrt( x*x + y*y + z*z );
      //cout << "dist : " << dist << endl;
      
      //Check if scatterers is far enough from others
      bool too_close = false;
      int i = 0;
      while( i < (int) array.size() && !too_close ){
	
	dist = sqrt( (x-array[i].x)*(x-array[i].x) + 
		     (y-array[i].y)*(y-array[i].y) + 
		     (z-array[i].z)*(z-array[i].z) );
	
	//cout << "(" << array[i].x << "," << array[i].y << "," << array[i].z << ")" << endl;
	//cout << "dist2 : " << dist << endl;
	
	if( dist < (2*r+d) )
	  too_close = true;
	
	i++;
      }
      
      if( !too_close ){
	//cout << "(" << x << "," << y << "," << z << ")" << endl;
	Pvec center(x, y, z, Pvec::Cartesian);
	array.push_back(c + center);
      }	
      
      
      n++;
  }
  
  return array;
  
  }




std::vector<Pvec> SphericalPeriodic(double R, Pvec center, double r, double p){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  //p: period length
  
  std::vector<Pvec> array; // Array of center locations

  int M = ceil( R / p );
  for( int i = -M; i <= M; i++ ){
    for( int j = -M; j <= M; j++ ){
      for( int k = -M; k <= M; k++ ){
	std::vector<double> Z(3);
	Z[0] = i*p;
	Z[1] = j*p;
	Z[2] = k*p;
	

	double rad = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );
	if( rad <= R ){
	  Pvec loc(Z[0] + center.x, Z[1] + center.y, Z[2] + center.z, Pvec::Cartesian);
	  array.push_back(loc);
	}
	
      }
    }
  }


  cout << "size : " << array.size() << endl;
  return array;  
}

//( Used only by ScatInit() )
// Returns a std::vector of location of sphere randomly located around the origin within
// a cubic enclosure
std::vector<Pvec> RandCubic(double Dx, double Dy, double Dz, double r, double d, Pvec c, int N){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere

  std::vector<Pvec> array; // Array of center locations

  while( true ){
    array.resize(0);
    
    //srand48(0);
    srand48(time(0));
    double dist = 0;
    double x,y,z;
    int maxN = 10*N; //Maximum number of iterations
    
    int n = 0;
    while( n < maxN  &&  (int) array.size() < N )
      {
	x = c.x + (2.*drand48() - 1.)*Dx;
	y = c.y + (2.*drand48() - 1.)*Dy;
	z = c.z + (2.*drand48() - 1.)*Dz;
	
	//Check sphere is far enough from others
	bool too_close = false;
	int i = 0;
	while( i < (int) array.size() ){
	  
	  dist = sqrt( (x-array[i].x)*(x-array[i].x) + 
		       (y-array[i].y)*(y-array[i].y) + 
		       (z-array[i].z)*(z-array[i].z) );
	  
	  if( dist < (2*r+d) ){
	    too_close = true;
	    break;
	  }
	  
	  //cout << "(" << x << "," << y << "," << z << ") : dist : " << dist << " : " << too_close << endl;
	  
	  i++;
	}
	
	if( !too_close ){
	  //cout << "(" << x << "," << y << "," << z << ")" << endl;
	  Pvec vec(x, y, z, Pvec::Cartesian);
	  Pvec loc = vec;
	  array.push_back(loc);
	}	
	
	n++;
      }

    if( (int) array.size() == (int) N ) { break; }
    //d *= 0.9;

    //cout << "d : " << d << endl;
  }
  
  return array;
  
}



// Samples uniform distribution within cube of side 2*R centered at the origin
std::vector<double> Uniform_Cube(double R)
{
  //srand48(0);
  //srand48(time(NULL));

  std::vector<double> x(3);
  for( int i = 0; i < 3; i++ )
    x[i] = 2*R*drand48() - R;
    
  return x;
}


// Samples uniform distribution within ball of radius R
std::vector<double> Uniform_Ball(double R)
{
  std::vector<double> x = Uniform_Cube(R);
  double d = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
  while( d > R )
    {
      x = Uniform_Cube(R);
      d = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
    }

  return x;
}


// Samples distribution within ball of radius (R+delta) which
// is uniform for [0,R) and decays smoothly (bump) in [R,R+delta)
double bumpF( double R, double delta, double x )
{
  if( std::abs(x) < R )
    return 1.;
  else if( std::abs(x) >= R && std::abs(x) <= (R+delta) )
    return exp(1.) * exp( -1. / (1 - pow((x-R)/delta, 2.)) );
  else
    return 0.;
}


std::vector<double> Bump_Ball(double R, double delta)
{
  double V = 4*PI*pow(R,3.)/3. + 4*PI*exp(1)*delta*( 0.2219969081*pow(R,2.) + 2.*0.0742477534*delta*R + 0.03510073838*pow(delta,2.) );
  //cout << "V : " << V << endl;

  std::vector<double> Z = Uniform_Ball(R+delta);
  double d = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );
  double U = drand48();
  while( U > (bumpF(R, delta, d)/V) )
    {
      Z = Uniform_Ball(R);
      d = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );
      U = drand48();
    }

  return Z;
}



// Cylinder of length L with axis parallel to the z-axis and radius R.
std::vector<Pvec> RandCyl(double L, double R, double r, double d, Pvec c, double N){

  srand48(0);

  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  int maxN = 100*N; //Maximum number of iterations
  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {

      std::vector<double> Z = Uniform_Cube(L);
      double D = sqrt( Z[0]*Z[0] + Z[1]*Z[1] );
      //cout << Z[1] << "," << Z[2] << "," << Z[3] << ";" << endl;
      
      if( D < R ){
	//Check sphere is far enough from others
	bool too_close = false;
	double dist;
	int i = 0;
	while( i < (int) array.size() && !too_close ){
	  
	  dist = sqrt( (Z[0]-array[i].x)*(Z[0]-array[i].x) + 
		       (Z[1]-array[i].y)*(Z[1]-array[i].y) + 
		       (Z[2]-array[i].z)*(Z[2]-array[i].z) );
	  
	  if( dist < (2*r+d) )
	    too_close = true;
	  
	  i++;
	}
	
	if( !too_close ){
	  Pvec center(Z[0], Z[1], Z[2], Pvec::Cartesian);
	  array.push_back(c + center);
	}	
      }
      
      n++;
    }
  
  return array;  
}



std::vector<Pvec> RandBump(double R, double delta, double r, double d, Pvec c, double N){

  srand48(8);

  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  int maxN = 120*N; //Maximum number of iterations
  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {

      std::vector<double> Z = Bump_Ball(R, delta);
      
      //Check sphere is far enough from others
      bool too_close = false;
      double dist;
      int i = 0;
      while( i < (int) array.size() && !too_close ){
	
	dist = sqrt( (Z[0]-(array[i].x-c.x))*(Z[0]-(array[i].x-c.x)) + 
		     (Z[1]-(array[i].y-c.y))*(Z[1]-(array[i].y-c.y)) + 
		     (Z[2]-(array[i].z-c.z))*(Z[2]-(array[i].z-c.z)) );
	

	if( dist < (2*r+d) )
	  too_close = true;
	

	i++;
      }
      
      if( !too_close ){
	Pvec center(Z[0], Z[1], Z[2], Pvec::Cartesian);
	array.push_back(c + center);
      }	
      
      n++;
    }
  
  cout << "NScat : " << array.size() << endl;
  return array;  
}






std::vector<Pvec> ErgodCavity2(double R, Pvec center, double r, double d, double N){

  srand48(8);
  //srand48(time(0));

  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere


  std::vector<Pvec> array; // Array of center locations

  double vol = 4*PI*pow(R,3.)/(3.*N);
  double side = pow(vol, 1./3.); 
  double sigma = std::max(0., side - 2*(r+d/2.));
  cout << "sigma  : " << sigma << endl;

  cout << " side : " << side << endl;
  
  // Periodic grid
  int M = ceil( R / side ) + 2;
  //int N = ceil( t / p ) + 2;
  for( int i = -M; i <= M; i++ ){
    for( int j = -M; j <= M; j++ ){
      for( int k = -M; k <= M; k++ ){

	// Center of cell
	std::vector<double> Z(3);
	Z[0] = (i+1./2.)*side;
	Z[1] = (j+1./2.)*side;
	Z[2] = (k+1./2.)*side;

	// Random perturbation
	Z[0] = (sigma*drand48() - sigma/2.) + Z[0];
	Z[1] = (sigma*drand48() - sigma/2.) + Z[1];
	Z[2] = (sigma*drand48() - sigma/2.) + Z[2];
	
	double rad = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );
	if( (Z[1]+center.y) < 0 && rad <= R ){
	  Pvec loc(Z[0] + center.x, Z[1] + center.y, Z[2] + center.z, Pvec::Cartesian);
	  array.push_back(loc);

	  cout << loc.x << "," << loc.y << "," << loc.z << ";";
	}
	
      }
    }
  }


  double MIN = 1e10;
  for( int i = 0; i < (int) array.size(); i++ ){
    for( int j = 0; j < (int) array.size(); j++ ){
      if( i != j )
	MIN = std::min(MIN, sqrt(pow(array[i].x-array[j].x, 2.) + pow(array[i].y-array[j].y, 2.) + pow(array[i].z-array[j].z, 2.)) );
    }
  }

  cout << "size : " << array.size() << endl;
  cout << "minimum distance : " << MIN << endl;
  return array;  
}


std::vector<Pvec> ErgodCavityPeriodic(double R, double t, Pvec center, double r, double p){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  //p: period length
  
  std::vector<Pvec> array; // Array of center locations

  int M = ceil( R / p ) + 2;
  //int N = ceil( t / p ) + 2;
  for( int i = -M; i <= M; i++ ){
    for( int j = -M; j <= M; j++ ){
      for( int k = -M; k <= M; k++ ){
	std::vector<double> Z(3);
	Z[0] = i*p;
	Z[1] = j*p;
	Z[2] = k*p;
	
	double rad = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );
	if( (Z[1]+center.y) < 0 && rad <= R ){ //&& std::abs(Z[2]) <= t){
	  Pvec loc(Z[0] + center.x, Z[1] + center.y, Z[2] + center.z, Pvec::Cartesian);
	  array.push_back(loc);
	}
	
      }
    }
  }


  cout << "size : " << array.size() << endl;
  return array;  
}

std::vector<Pvec> ErgodCavityPeriodicRandom(double R, double t, Pvec center, double r, double d, double p){
  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  //p: period length

  double tau = 0.1*p;
  srand48(3);
  
  std::vector<Pvec> array; // Array of center locations

  int M = ceil( R / p ) + 2;
  int N = ceil( t / p ) + 2;
  for( int i = -M; i <= M; i++ ){
    for( int j = -M; j <= M; j++ ){
      for( int k = -N; k <= N; k++ ){
	std::vector<double> Z(3);
	Z[0] = i*p;
	Z[1] = j*p;
	Z[2] = k*p;
	
	double rad = sqrt( Z[0]*Z[0] + Z[1]*Z[1] );
	if( (Z[1]+center.y) < 0 && rad <= R && std::abs(Z[2]) <= t){
	  Pvec loc(Z[0] + tau*drand48() + center.x, 
		   Z[1] + tau*drand48() + center.y, 
		   Z[2] + tau*drand48() + center.z, Pvec::Cartesian);
	  array.push_back(loc);
	}
	
      }
    }
  }


  cout << "size : " << array.size() << endl;
  return array;  
}





std::vector<Pvec> ErgodShell(double R, Pvec center, double r, double d, double N){

  srand48(3);

  //r is the radius of the spheres
  //d is the std::minimum distance allowable between any sphere
  
  std::vector<Pvec> array; // Array of center locations
  int maxN = 100*N; //Maximum number of iterations
  int n = 0;
  while( n < maxN  &&  (int) array.size() < N )
    {

      // Create random point inside sphere
      std::vector<double> Z(3);
      Z[0] = (2*R*drand48() - R);
      Z[1] = (2*R*drand48() - R);
      Z[2] = (2*R*drand48() - R);

      //cout << Z[0] << "," << Z[1] << "," << Z[2] << endl;


      //Check sphere is far enough from others
      bool too_close = false;
      double dist;
      int i = 0;
      // Creates edge
      double rad = sqrt( Z[0]*Z[0] + Z[1]*Z[1] + Z[2]*Z[2] );

      //cout << Z[0]-center.x  << endl;

      if( rad <= R){
	//&& std::abs(Z[0]) > 1e-10
	// && std::abs(Z[1]) > 1e-10
	// && std::abs(Z[2]) > 1e-10){
	
	// Project on shell
	//Z[0] *= R/rad; 
	Z[0] += center.x; 
	//Z[1] *= R/rad; 
	Z[1] += center.y; 
	//Z[2] *= R/rad; 
	Z[2] += center.z; 

	// If y-component positive, project on x-z plane
	if( Z[1] > 0 )
	  Z[1] = 0;

	// Check if not too close to another sphere
	while( i < (int) array.size() && !too_close ){
	  
	  dist = sqrt( (Z[0]-array[i].x)*(Z[0]-array[i].x) + 
		       (Z[1]-array[i].y)*(Z[1]-array[i].y) + 
		       (Z[2]-array[i].z)*(Z[2]-array[i].z) );
	  
	  
	  if( dist < (2*r+d) )
	    too_close = true;
	  
	  i++;
	}
	
	if( !too_close ){
	  Pvec loc(Z[0], Z[1], Z[2], Pvec::Cartesian);
	  array.push_back(loc);
	}	
      }

      n++;
    }
  
  cout << "Total number of scatterers : " << array.size() << endl;
  return array;  
}


#endif
