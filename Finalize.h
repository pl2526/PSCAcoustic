/*
 *  Finalize.h
 *  PSCAcoustic
 *
 *  Post-processing for writing solution to file in 
 *  standard form.
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


#ifndef FINALIZE_H
#define FINALIZE_H

#include "Constants.h"
#include "General.h"
#include "Coordinates.h"
#include "Scatterer.h"
#include "IncWave.h"

// Write to file information pertaining to the problem solved
void Write_Info(int Proc_Idx, std::string& type,double residual, double relative_residual, double cond, int Niter, std::string filename, complex k_out = K_OUT, complex k_in = K, int N = NScat);

// Write to file information about each source
void Write_Source( int Proc_Idx, IncWave* IW, IncWave::Type type, std::string filename);
void Write_Source( int Proc_Idx, std::vector<IncWave*>& IW, IncWave::Type type, std::string filename);

// Write to file location of each scatterer
void Write_Location( int Proc_Idx, std::vector<Scatterer>& ScL, std::string filename);

// Write to file solution in VSWF basis
void Write_Solution( int Proc_Idx, Indexing& index, std::vector< std::vector<complex> > u, std::string filename);


#endif
