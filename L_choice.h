/*
 *  L_choice.h
 *  PSCAcoustic
 *
 *  Routines used to established the number of coefficients necessary for each scatterers
 *  and at each level in order to achieve a certain accuracy.
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

#ifndef L_CHOICE_H
#define L_CHOICE_H

#include "General.h"
#include "Coordinates.h"
#include "Scatterer.h"
#include "./FMMPS/TransferUtils.h"

int L_max_scat(double boxSize, double d_min, double eps, complex k, complex k_out, int N = 40);     // Parameter (Largest degree os SWF) for scatterers
int L_max_level(double boxSize, double eps, double radius, complex k, complex k_out, int N = 50 );       // Parameter (Largest degree os SWF) for tree levels


#endif
