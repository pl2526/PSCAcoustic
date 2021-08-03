/*
 *  Quadrature.h
 *  PSCAcoustic
 *
 *  Quadrature for High-frequency portion of code
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


#ifndef QUADRATURE_FMM_H
#define QUADRATURE_FMM_H

#include "General.h"
#include "Vec3.h"
#include "HelmholtzUtil.h"
#include "../Scatterer.h"



namespace FMM {

  // A quadrature point 
  struct QuadraturePoint
  {
    int index;
    double theta, phi, w;
    Vec3 s;

    QuadraturePoint( double theta_ = 0, double phi_ = 0, double w_ = 0 );
    ~QuadraturePoint() {}
  };
  

  
  // A row of quadrature points on the sphere
  struct QuadratureRow
  {
    int index;
    double theta;
    std::vector<QuadraturePoint> point;
    
    // Constructor
    // A row of N_phi quadrature points (of total weight w0) at constant theta
    QuadratureRow( double theta_ = 0, int N_phi = 0, double w0 = 0 );

    // Destructor
    ~QuadratureRow() {}

    // Accessors    
    inline int size() { return point.size(); }
    inline int numPoints() { return size(); }
    inline double getTheta() { return theta; }
    inline QuadraturePoint& getPoint( int k ) { return point[k]; }
    
    // Output
    friend ostream& operator<<(ostream& os, const QuadratureRow& quad){ return os; }
  };
  
  

  class QuadratureHelm
  {
  public:
    // TODO : Should it be complex?
    const double kappa;
    
    std::vector<QuadratureRow> row;         // List of quadrature rows (const theta)
    std::vector<QuadraturePoint*> point;    // Enumerated quadrature points
    
    int L;                             // Gegenbauer truncation
    int max_N_phi;
    
    // Index maps for rotations and reflections of the quadrature
    std::vector< std::vector<int> > indexMap;
    
    QuadratureHelm( HelmKernel& K );

    // Destructor
    ~QuadratureHelm(){}
    
    // Accessors
    inline int size() const { return point.size(); }
    inline int numPoints() const { return size(); }
    inline int numRows() const { return row.size(); }
    inline int maxRowSize() const { return max_N_phi; }
    inline int getTruncation() const { return L; }
    inline QuadratureRow& getRow( int k ) { return row[k]; }
    inline QuadraturePoint& getPoint( int k ) { return *(point[k]); }
    
    // Output
    friend ostream& operator<<(ostream& os, const QuadratureHelm& q)
   {
      os << "(" << q.getTruncation() 
	 << "," << q.numRows() 
	 << "," << q.size() << ")" << endl;
      return os;
    }
    

  protected:
    
    inline void makePointerAccess()
    {
      // Quadrature is done being constructed
      // Construct the enumerated point std::vector
      int rowNumber = 0, pointNumber = 0;
      max_N_phi = 0;
      for( int k = 0; k < numRows(); ++k ) {
	// Assign each a row number
	row[k].index = rowNumber++;
	// Find the largest row
	max_N_phi = std::max( max_N_phi, row[k].size() );
	
	for( int p = 0; p < row[k].numPoints(); ++p ) {
	  // Assign each a point number
	  row[k].point[p].index = pointNumber++;
	  // Insert point
	  point.push_back( &(row[k].point[p]) );
	}
      }
    }

  };
  
  

  // A uniform quadrature in phi and theta
  class Quadrature_Uniform : public QuadratureHelm
  {
    complex k;
    complex k_out;

  public:

    // Constructor
    Quadrature_Uniform(HelmKernel& K, double boxSize, double eps, complex k, complex k_out);
    
    // Destructor
~Quadrature_Uniform() {}

//----- Auxilliary functions -----//
double g( int n, double C, double k_out, double r);
double h( int n, double C, double k_out, double r, double theta);

    //----- Functions providing parameters values -----//
int get_N_theta( int L, double k_out, double r0, double r, double boxSize, double eps );
int get_N_phi( std::vector< std::vector<complex> >& TL, 
		 int n,  int N_phiT, double maxL, double theta, 
		 int N_theta, double k_out, double r, double boxSize, double eps);
  
  };
  
  
  typedef Quadrature_Uniform Quadrature;
  
}

#endif
