/*
 *  Retorpolator.h
 *  PSCAcoustic
 *
 *  Fourier interpolation and anterpolation between levels
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


#ifndef INTERPOLATOR_FMM_H
#define INTERPOLATOR_FMM_H

#include "Quadrature.h"
#include "SphericalHarmonicsTransform.h"
#include "General.h"

namespace FMM {
  
  
  // A reterpolator between two quadratures using FFT interpolation
  // Optimized to reduce the number of FFTs
  // Generalized to handle any quadrature
  // only used in the HF region
  // 
  // ***Always goes from 1 to 2***
  class FFTInterp_FFT2
  {
    
  public:
    Quadrature* q1;
    Quadrature* q2;
    
  private:
    // q1 Data
    int N_rows1, N_theta1, N_phi1;
    
    // q2 Data
    int N_rows2, N_theta2, N_phi2;
    
    // Intermediate Data
    int N_rows, N_phi, N_phi0;
    
    // Intermediate Storage  N_rows x N_phi
    std::vector<complex> T;
    // Intermediate Storage  N_theta
    std::vector<complex> t;
    
    // First stage phi (row) FFTs
    std::vector<fftw_plan> rFFT1;
    
    // Second stage theta (column) FFTs
    fftw_plan tFFT;
    fftw_plan tIFFT;
    
    // Third stage phi (row) FFTs
    std::vector<fftw_plan> rIFFT2;
    

  public:
    
    // Constructor
    // Interpolation from q1 to q2
    FFTInterp_FFT2( Quadrature* q1_, Quadrature* q2_);
    
    // Destructor
    ~FFTInterp_FFT2();
    
    
    // HF Interpolation from A(1) to B(2)
    void apply( std::vector<complex>& vec_A, Quadrature* quad_A,
		std::vector<complex>& vec_B, Quadrature* quad_B);
    
  private:
    // Disable Copy and Assignment
    FFTInterp_FFT2( const FFTInterp_FFT2& S ) {}
    void operator=( const FFTInterp_FFT2& S ) {}
  };
  
  typedef class FFTInterp_FFT2 Interpolator;
  typedef class FFTInterp_FFT2 Anterpolator;
  

}


#endif
