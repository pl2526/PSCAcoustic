/*
 *  Coordinates.h
 *  PSCAcoustic
 *
 *  Generalized 3D std::vectors containing coordinates in both Spherical and
 *  Cartesian reference frames.
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
