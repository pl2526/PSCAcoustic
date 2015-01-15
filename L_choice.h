/*
 *  L_choice.h
 *  PSCAcoustic
 *
 *  Routines used to established the number of coefficients necessary for each scatterers
 *  and at each level in order to achieve a certain accuracy.
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

#ifndef L_CHOICE_H
#define L_CHOICE_H

#include "General.h"
#include "Coordinates.h"
#include "Scatterer.h"
#include "./FMMPS/TransferUtils.h"

int L_max_scat(double boxSize, double d_min, double eps, complex k, complex k_out, int N = 40);     // Parameter (Largest degree os SWF) for scatterers
int L_max_level(double boxSize, double eps, double radius, complex k, complex k_out, int N = 50 );       // Parameter (Largest degree os SWF) for tree levels


#endif
