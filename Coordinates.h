/*
 *  Coordinates.h
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


#ifndef COORDINATES_H
#define COORDINATES_H

#include "Constants.h"

//Prototypes
void Cart2Sf( double x, double y, double z, double& r, double& theta, double& phi );
void Sf2Cart( double r, double theta, double phi, double& x, double& y, double& z );

//Pvec : position std::vector
//Structure representing a point in 3D
//contains its coordinate in both cartesian and spherical coordinates
struct Pvec{

  enum Type{Cartesian, Spherical};

  double x,y,z; //Cartesian
  double r,theta,phi;  //Spherical

  //Empty constructors
  Pvec(){}
  Pvec(double p1, double p2, double p3, Type T );
    
  //Empty destructor
  ~Pvec(){}

  //Overload +/- operations on Geeneralized Vectors
  friend Pvec operator + (Pvec vec1, Pvec vec2);
  friend Pvec operator - (Pvec vec1, Pvec vec2);
  friend Pvec operator - (Pvec vec1);
  friend Pvec operator * (double a, Pvec vec1);
  //Overload == operations on Generalized Vectors
  friend bool operator == (Pvec vec1, Pvec vec2);
  
};




#endif
