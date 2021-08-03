/*
 *  Vec3.h
 *  PSCAcoustic
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

#ifndef VEC3_FMM_H
#define VEC3_FMM_H

namespace FMM {

  struct Vec3
  {
    double x, y, z;
    
    // Constructor
  Vec3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
    // Copy Constructor
  Vec3(const Vec3& v) : x(v.x), y(v.y), z(v.z) {}
    
    /**********************/
    // Operator Overloads //
    /**********************/
    
    inline bool operator==(const Vec3& v) const {
      return (x == v.x && y == v.y && z == v.z);
    }
    
    inline Vec3 operator+(const Vec3& v) const { 
      return Vec3(x + v.x, y + v.y, z + v.z);
    }
    
    inline Vec3 operator+(double S) const {
      return Vec3(x + S, y + S, z + S);
    }
    
    inline void operator+=(const Vec3& v) {
      x += v.x; y += v.y; z += v.z;
    }
    
    inline void operator+=(double S) {
      x += S; y += S; z += S;
    }
    
    inline Vec3 operator-() const {
      return Vec3(-x, -y, -z);
    }
    
    inline Vec3 operator-(const Vec3& v) const {
      return Vec3(x - v.x,  y - v.y,  z - v.z);
    }
    
    inline Vec3 operator-(double S) const {
      return Vec3(x - S, y - S, z - S);
    }
    
    inline void operator-=(const Vec3& v) {
      x -= v.x; y -= v.y; z -= v.z;
    }
    
    inline void operator-=(double S) {
      x -= S; y -= S; z -= S;
    }
    
    inline Vec3 operator/(const Vec3& v) const {
      return Vec3(x / v.x,  y / v.y,  z / v.z);
    }
    
    inline Vec3 operator/(double S) const {
      return Vec3(x / S, y / S, z / S);
    }
    
    inline void operator/=(const Vec3& v) {
      x /= v.x; y /= v.y; z /= v.z;
    }
    
    inline void operator/=(double S) {
      x /= S; y /= S; z /= S;
    }
    
    inline Vec3 operator*(const Vec3& v) const {
      return Vec3(x * v.x,  y * v.y,  z * v.z);
    }
    
    inline Vec3 operator*(double S) const {
      return Vec3(x * S,  y * S,  z * S);
    }
    
    inline void operator*=(const Vec3& v) {
      x *= v.x; y *= v.y; z *= v.z;
    }
    
    inline void operator*=(double S) {
      x *= S; y *= S; z *= S;
    }
    
    /*************/
    // Functions //
    /*************/
    
    inline double dot(const Vec3 &v) const {
      return x*v.x + y*v.y + z*v.z;
    }
    
    inline Vec3 cross(const Vec3 &v) const {
      return Vec3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    
    inline double magSq() const {
      return x*x + y*y + z*z;
    }
    
    inline double mag() const {
      return sqrt( this->magSq() );
    }
    
    inline double dist(const Vec3& v) const {
      return (v - *this).mag();
    }
    
    inline void normalize() {
      double mag = this->mag();
      if(mag != 0) (*this) *= 1.0/mag;
    }
    
    /**********/
    // Output //
    /**********/
    
    friend ostream& operator<<(ostream& os, const Vec3& v) {
      return (os << "(" << v.x << "," << v.y << "," << v.z << ")");
    }
  };
  
}

#endif
