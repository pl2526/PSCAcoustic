/*
 *  Finalize.cpp
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

#ifndef FINALIZE_CPP
#define FINALIZE_CPP

#include "Finalize.h"


void Write_Info(int Proc_Idx, std::string& type,double residual, double relative_residual, double cond, int Niter, std::string filename, complex k_out, complex k_in, int N){
      std::stringstream nProc_str;
      nProc_str << Proc_Idx;
      filename.append(nProc_str.str());
      filename.append("_Inf.csv");
      char *Filename = (char*) filename.c_str();
      ofstream Infofile(Filename, ios::out);
      
      // Index
      Infofile << " Simulation index :," << Proc_Idx << endl << endl;

      // Linear solver parameters
      Infofile << " Relative error tolerance :," << RTOL << endl;
      Infofile << " Absolute error tolerance :," << ATOL << endl;
      Infofile << " Divergence tolerance :," << DTOL << endl;
      Infofile << " Maximum number of iterations :," << MAXITS << endl << endl;
      
      // Background properties
      Infofile << " Background medium wave numver (K_OUT) :," << std::real(k_out) << "," << std::imag(k_out) << endl;
      Infofile << " Background medium density (RHO_OUT) :," << RHO_OUT << endl << endl;
      
      // Scatterer properties
      Infofile << " Scatterer wave numver (K) :," << std::real(k_in) << "," << std::imag(k_in) << endl;
      Infofile << " Scatterer density (RHO) :," << RHO << endl;
      Infofile << " Scatterer radius (RADIUS) :' ," << RADIUS << endl;
      Infofile << " Total number of scatterers (NScat) :," << N << endl;
      Infofile << " Total number of unknowns :," << NScat*(LMAX+1)*(LMAX+1)  << endl;
      Infofile << " Number of unknowns per scatterer :," << (LMAX+1)*(LMAX+1) << endl << endl;

      // Fast algorithm parameters
      Infofile << " Maximum degree of SWF (Cartesian) expansion (LMAX) :," << LMAX << endl;
      Infofile << " Number of levels in FMM trees (nLevels) :," << nLevels << endl;
      Infofile << " FMM accuracy (eps) :," << EPS << endl;
      Infofile << " Interface between high and low frequency (CUTOFF_SIZE) :," << CUTOFF_SIZE << endl << endl;
	  
      // Solution characteristics
      Infofile << " L-S residual :," << std::setprecision(15) << residual << endl;
      Infofile << " L-S relative residual :," << std::setprecision(15) << relative_residual << endl;
      Infofile << " Approximate condition number :," << std::setprecision(15) << cond << endl;
      Infofile << " Number of iterations :," << Niter << endl << endl;

      // TODO: necessary?
      // flags
      Infofile << " CT_PRECOMPUTE :," << CT_PRECOMPUTE << endl;
      Infofile << " FAST_TRANSFER :," << FAST_TRANSFER << endl << endl;

      Infofile.close();

}


void Write_Source( int Proc_Idx,IncWave* IW, IncWave::Type type, std::string filename){
  std::vector<IncWave*> IW_vec(1); IW_vec[0] = IW;
  Write_Source(Proc_Idx, IW_vec, type, filename);
}
  
void Write_Source( int Proc_Idx, std::vector<IncWave*>& IW, IncWave::Type type, std::string filename){
  
  std::stringstream nProc_str;
  nProc_str << Proc_Idx;
  filename.append(nProc_str.str());
  filename.append("_Src.csv");
  char *Filename = (char*) filename.c_str();
  ofstream Srcfile(Filename, ios::out);
  


  if( type == IncWave::PW ){
    

   // Sorce properties
    /* Infofile << " Source type :," << type << endl;
      Infofile << " Direction(if plane wave)/ location(if point source) : ', " << d.x << "," << d.y << "," << d.z << endl;
      Infofile << " Longitudinal wave component :, " << long_comp << endl;
      Infofile << " Shear wave component :, " << perp_comp.x << "," << perp_comp.y << "," << perp_comp.z << endl << endl;
      */
    
  } else if( type == IncWave::Pt ){

    Srcfile << "Source type:, Pt" << endl;
    Srcfile << "Idx , amplitude, x-location, y-location, z-location" << endl;

    for( int n = 0; n < (int) IW.size(); n++ ){
      SphericalWave* IW_ = static_cast<SphericalWave*>(IW[n]);
      Srcfile << n << " , " << std::real(IW_->getAmplitude()) <<  " , " << IW_->getLocation().x << 
	" , " << IW_->getLocation().y << " , " << IW_->getLocation().z << endl;
    }
  }

}

void Write_Location( int Proc_Idx, std::vector<Scatterer>& ScL, std::string filename){
  
  //std::string filename(LOC_PREFIX);
  std::stringstream nProc_str;
  nProc_str << Proc_Idx;
  filename.append(nProc_str.str());
  filename.append("_Loc.csv");
  char *Filename = (char*) filename.c_str();
  ofstream Locfile(Filename, ios::out);
  
  for( int n = 0; n < (int) ScL.size(); n++ ){
    Pvec p = ScL[n].location;
    Locfile << n << ","  << p.x << "," << p.y << "," << p.z << endl;
  }

}

// Write components of solution in VSWF basis
void Write_Solution( int Proc_Idx, Indexing& index, std::vector< std::vector<complex> > u, std::string filename){

  //std::string filename(SOL_PREFIX);
  std::stringstream nProc_str;
  nProc_str << Proc_Idx;
  filename.append(nProc_str.str());
  filename.append("_Sol.csv");
  char *Filename = (char*) filename.c_str();
  ofstream Solfile(Filename, ios::out);
  
  for( int n = 0; n < (int) u.size(); n++ )
    for( int i = 0; i < (int) index.size(); i++ )
      Solfile << n << ","  << index(i,0) << "," << index(i,1) << "," << std::real(u[n][i]) << "," << std::imag(u[n][i]) << endl;


}








#endif
