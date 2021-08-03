/*
 *  Coordinates.cpp
 *  PSCAcoustic
 *
 *  Generalized 3D std::vectors containing coordinates in both Spherical and
 *  Cartesian reference frames.
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


#ifndef COORDINATES_CPP
#define COORDINATES_CPP

#include <vector>
#include <cmath>
#include <complex>
#include <assert.h>

#include "Coordinates.h"
//#include "Constants.h"

// Constructor
Pvec::Pvec(double p1, double p2, double p3, Type T ){
  
  if ( T == Cartesian ){
    x = p1;
    y = p2;
    z = p3;
    Cart2Sf(x, y, z, r, theta, phi); //Conversion
    
  } else if ( T == Spherical ){
    
    assert( p1 >= 0 );
    r = p1;
    theta = p2;
    phi = p3;
    Sf2Cart(r, theta, phi, x, y, z); //Conversion
    
  }   
  
}


//Overload +/- operations on Generalized Vectors
Pvec operator + (Pvec vec1, Pvec vec2){
  double vec0_x = vec1.x + vec2.x;
  double vec0_y = vec1.y + vec2.y;
  double vec0_z = vec1.z + vec2.z;
  
  return Pvec(vec0_x, vec0_y, vec0_z, Pvec::Cartesian);
}

Pvec operator - (Pvec vec1, Pvec vec2){
  double vec0_x = vec1.x - vec2.x;
  double vec0_y = vec1.y - vec2.y;
  double vec0_z = vec1.z - vec2.z;
  
  return Pvec(vec0_x, vec0_y, vec0_z, Pvec::Cartesian);
}

Pvec operator - (Pvec vec1){
  double vec0_x = -vec1.x;
  double vec0_y = -vec1.y;
  double vec0_z = -vec1.z;
  
  return Pvec(vec0_x, vec0_y, vec0_z, Pvec::Cartesian);
}

Pvec operator * (double a, Pvec vec1){
  double vec0_x = a*vec1.x;
  double vec0_y = a*vec1.y;
  double vec0_z = a*vec1.z;
  
  return Pvec(vec0_x, vec0_y, vec0_z, Pvec::Cartesian);
}

//Overload == operations on Geeneralized Vectors
bool operator == (Pvec vec1, Pvec vec2){
  if( (vec1.x == vec2.x) && (vec1.y == vec2.y) && (vec1.z == vec2.z) )
    return true;
  else
    return false;
}





//***Auxiliary functions***

//Goes from cartesian coordinates to spherical coordinates
void Cart2Sf( double x, double y, double z, double& r, double& theta, double& phi ){
  
  r = sqrt( x*x + y*y + z*z );
  if ( r < 1e-15 ) {
    theta = 0; phi = 0;
    return;
  } else {
    theta = acos(z/r);
    
    //Phi is undefined if location is at a pole
    if ( std::abs(theta-PI) < 1e-15 || std::abs(theta) < 1e-15 ){
      phi = 0;
    } else {
      phi = atan( y/x );
      if( x < 0 ) { phi += PI; }
    }
  }
  
  return;
}


//Goes from spherical coordinates to cartesian coordinates 
void Sf2Cart( double r, double theta, double phi, double& x, double& y, double& z ){
  
  assert( theta >= 0 );
  assert( theta <= PI );
  
  x = r * sin(theta) * cos(phi);
  y = r * sin(theta) * sin(phi);
  z = r * cos(theta);
  
  return;
}




#endif
