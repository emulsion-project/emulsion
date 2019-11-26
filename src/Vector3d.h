#ifndef VECTOR3D_H_INCLUDED
#define VECTOR3D_H_INCLUDED

#include <ostream>
#include "Complex.h"

/*********************************
 *     Class declaration         *
 *********************************/

class Vector3d {
public:
	Vector3d();
	Vector3d(const Complex &x);
	Vector3d(const Vector3d &v);
	~Vector3d();
	Vector3d& operator=(const Vector3d &v);
	Complex& operator()(const int &i);
	Complex& operator()(const int &i) const;
	int n;
	Complex *t;
};

/*********************************
 *     Associated operators      *
 *********************************/

Vector3d operator+(const Vector3d &v1, const Vector3d &v2);
Vector3d operator-(const Vector3d &v1, const Vector3d &v2);
Vector3d operator*(const Vector3d &v1, const Complex &z);
Vector3d operator*(const Complex &z, const Vector3d &v1);
std::ostream& operator<<(std::ostream &out, const Vector3d &v);

/*********************************
 *     Associated functions      *
 *********************************/

Complex dotProduct(const Vector3d &v1, const Vector3d &v2);
Vector3d crossProduct(const Vector3d &v1, const Vector3d &v2);
Vector3d toSphericalCoordinates(const Vector3d &v, const double &t,
		const double &p);
Vector3d toCartesianCoordinates(const Vector3d &v, const double &t,
		const double &p);
Vector3d toSphericalCoordinates(const Vector3d &v, const Complex &t,
		const double &p);
Vector3d toCartesianCoordinates(const Vector3d &v, const Complex &t,
		const double &p);

#endif // VECTOR3D_H_INCLUDED
