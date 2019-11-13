#include "Vector3d.h"

/*********************************
 *     Class constructors        *
 *********************************/

Vector3d::Vector3d() : n(3), t(new Complex[3]) {}
Vector3d::Vector3d(const Complex &x) : n(3), t(new Complex[3]) { for (int i=0 ; i<n ; i++) t[i] = x; }
Vector3d::Vector3d(const Vector3d &v) : n(3), t(new Complex[3]) { for (int i=0 ; i<n ; i++) t[i] = v.t[i]; }

/*********************************
 *     Class destructor          *
 *********************************/

Vector3d::~Vector3d() { delete[] (t); }

/*********************************
 *     Class member operators    *
 *********************************/

Vector3d& Vector3d::operator=(const Vector3d &v) {
  if (this != &v) {
    for (int i=0 ; i<n ; i++) t[i] = v.t[i];
  }
  return *this;
}

Complex& Vector3d::operator()(const int &i) { return t[i]; }
Complex& Vector3d::operator()(const int &i) const { return t[i]; }

/*********************************
 *     Associated operators      *
 *********************************/

Vector3d operator+(const Vector3d &v1, const Vector3d &v2) {
  Vector3d v3;
  for (int i=0 ; i<3 ; i++) v3.t[i] = v1.t[i]+v2.t[i];
  return v3;
}

Vector3d operator-(const Vector3d &v1, const Vector3d &v2) {
  Vector3d v3;
  for (int i=0 ; i<3 ; i++) v3.t[i] = v1.t[i]-v2.t[i];
  return v3;
}

Vector3d operator*(const Vector3d &v1, const Complex &z) {
  Vector3d v2;
  for(int i=0 ; i<3 ; i++) v2.t[i] = v1.t[i] * z;
  return v2;
}

Vector3d operator*(const Complex &z, const Vector3d &v1) {
  Vector3d v2;
  for(int i=0 ; i<3 ; i++) v2.t[i] = v1.t[i] * z;
  return v2;
}

std::ostream &operator<<(std::ostream &out, const Vector3d &v) {
  for(int i=0 ; i<3 ; i++) out << v.t[i] << "\n";
  return out;
}

/*********************************
 *     Associated functions      *
 *********************************/

Complex dotProduct(const Vector3d &v1, const Vector3d &v2) {
  Complex z(0.);
  for (int i=0 ; i<3 ; i++) z += (v1.t[i])*conj(v2.t[i]);
  return z;
}

Vector3d crossProduct(const Vector3d &v1, const Vector3d &v2) {
  Vector3d v;
  v.t[0] = v1.t[1]*v2.t[2] - v1.t[2]*v2.t[1];
  v.t[1] = v1.t[2]*v2.t[0] - v1.t[0]*v2.t[2];
  v.t[2] = v1.t[0]*v2.t[1] - v1.t[1]*v2.t[0];
  return v;
}

Vector3d toCartesianCoordinates(const Vector3d &v, const double &t, const double &p) {
  double cost = cos(t), sint = sin(t), cosp = cos(p), sinp = sin(p);
  Vector3d v1(0.);
  v1.t[0] = sint*cosp*v.t[0] + cost*cosp*v.t[1] - sinp*v.t[2];
  v1.t[1] = sint*sinp*v.t[0] + cost*sinp*v.t[1] + cosp*v.t[2];
  v1.t[2] = cost*v.t[0] - sint*v.t[1];
  return v1;
}

Vector3d toSphericalCoordinates(const Vector3d &v, const double &t, const double &p) {
  double cost = cos(t), sint = sin(t), cosp = cos(p), sinp = sin(p);
  Vector3d v1(0.);
  v1.t[0] = sint*cosp*v.t[0] + sint*sinp*v.t[1] + cost*v.t[2];
  v1.t[1] = cost*cosp*v.t[0] + cost*sinp*v.t[1] - sint*v.t[2];
  v1.t[2] = -sinp*v.t[0] + cosp*v.t[1];
  return v1;
}

Vector3d toCartesianCoordinates(const Vector3d &v, const Complex &t, const double &p) {
  double cosp = cos(p), sinp = sin(p);
  Complex cost = cos(t), sint = sin(t);
  Vector3d v1(0.);
  v1.t[0] = sint*cosp*v.t[0] + cost*cosp*v.t[1] - sinp*v.t[2];
  v1.t[1] = sint*sinp*v.t[0] + cost*sinp*v.t[1] + cosp*v.t[2];
  v1.t[2] = cost*v.t[0] - sint*v.t[1];
  return v1;
}

Vector3d toSphericalCoordinates(const Vector3d &v, const Complex &t, const double &p) {
  double cosp = cos(p), sinp = sin(p);
  Complex cost = cos(t), sint = sin(t);
  Vector3d v1(0.);
  v1.t[0] = sint*cosp*v.t[0] + sint*sinp*v.t[1] + cost*v.t[2];
  v1.t[1] = cost*cosp*v.t[0] + cost*sinp*v.t[1] - sint*v.t[2];
  v1.t[2] = -sinp*v.t[0] + cosp*v.t[1];
  return v1;
}
