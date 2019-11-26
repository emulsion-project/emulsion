#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED

#include <ostream>
#include "Complex.h"

/*********************************
 *     Class declaration         *
 *********************************/

class Vector {
public:
	Vector();
	Vector(const int &n);
	Vector(const int &n, const Complex &x);
	Vector(const Vector &v);
	~Vector();
	void resize(const int &m);
	void assign(const int &m, const Complex &x);
	Vector& operator=(const Vector &v);
	int n;
	Complex *t;
};

/*********************************
 *     Associated operators      *
 *********************************/

Vector operator+(const Vector &v1, const Vector &v2);
Vector operator-(const Vector &v1, const Vector &v2);
Complex operator*(const Vector &v1, const Vector &v2);
Vector operator*(const Vector &v1, const Complex &z);
Vector operator*(const Complex &z, const Vector &v1);
std::ostream& operator<<(std::ostream &out, const Vector &v);

/*********************************
 *     Associated functions      *
 *********************************/

Complex max(const Vector &v);
Complex min(const Vector &v);
Complex sum(const Vector &v);

Vector poleAmplitude(const Vector &x, const Vector &y);
double norm(const Vector &v);
Vector fft(const Vector &v, const int &isign);

#endif // VECTOR_H_INCLUDED
