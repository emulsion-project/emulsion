#ifndef ROTATION_H_INCLUDED
#define ROTATION_H_INCLUDED

#include "WignerFunctions.h"
#include "Matrix.h"

class Rotation {
public:
	Rotation();
	Rotation(const int &n, const double &a, const double &b, const double &g);
	void setParameters(const int &n, const double &a, const double &b,
			const double &g);
	void calcRotationMatrices();
	Vector rotateVectorD(const Vector &v);
	Vector rotateVectorI(const Vector &v);
	Matrix getRd();
	Matrix getRi();

	int N, NM;
	double alpha, beta, gamma;
	Matrix Rd, Ri;
	WignerFunctions wigner;
};

#endif // ROTATION_H_INCLUDED
