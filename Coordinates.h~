/*
 *  Coordinates.h
 *  PSC
 *
 *  Created by Pierre-David Letourneau on 1/9/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 *  Generalized 3D std::vectors containing coordinates in both Spherical and
 *  Cartesian reference frames.
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
