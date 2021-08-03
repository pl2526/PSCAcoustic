/*
 *  Imaging.h
 *  PSCAcoustic
 *
 *  Produces image of wave field.
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


#ifndef IMAGING_H
#define IMAGING_H

#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "IncWave.h"
#include "Scatterer.h"

// u_in: RSWF expansion for field inside scatterers

void Imaging( std::vector<Scatterer>& ScL, IncWave* IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL, complex k_out = K_OUT, 
	      complex k_in = K);

void Imaging( std::vector<Scatterer>& ScL, std::vector<IncWave*>& IW, std::string filename,
	      std::vector< std::vector<complex> >& u, int idx, int M, double L,
	      std::vector< std::vector<complex> >* u_in = NULL, complex k_out = K_OUT, 
	      complex k_in = K);


#endif
