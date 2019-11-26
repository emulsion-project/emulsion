#ifndef TRANSLATION_H_INCLUDED
#define TRANSLATION_H_INCLUDED

#include "Rotation.h"
#include "ClebshGordan.h"
#include "LegendreFunctions.h"
#include "BesselFunctions.h"

class Translation {
public:
	Translation();
	Translation(const int &typ, const bool &dir, const int &n, const Complex kk,
			const double &xx, const double &yy, const double &zz);
	void setParameters(const int &typ, const bool &dir, const int &n,
			const Complex kk, const double &xx, const double &yy,
			const double &zz);
	void calcMatrix();
	void calcMatrixZ();
	Vector translate(const Vector &v, const int &direction);
	Vector translateT(const Vector &v, const int &direction);
	Matrix getMatrix(const int &direction);
	Matrix getMatrixT(const int &direction);

	ClebshGordan cg;
	SphericalBessel bes;
	LegendreFunctions leg;
	WignerFunctions wigner;
	int N, NM, type;
	bool direct;
	double x, y, z, r, t, p;
	Matrix Ad, Ai, Bd, Bi;
	Complex k;
};

#endif // TRANSLATION_H_INCLUDED
