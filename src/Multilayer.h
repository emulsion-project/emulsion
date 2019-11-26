#ifndef MULTILAYER_H_INCLUDED
#define MULTILAYER_H_INCLUDED

#include "Matrix.h"
#include "RealContainers.h"
#include <fstream>
#include "RefractiveIndex.h"
#include "Vector3d.h"

class Multilayer {
public:
	Multilayer();
	Multilayer(const string &fname, const double &lambd);
	void setParameters(const string &fname, const double &lambd);
	void setSourceExt(const Complex &e0, const Complex &beta,
			const double &alpha);
	Vector setSourceInt(const int &pol, const Complex &aP, const Complex &aM,
			const Complex &betaP, const double &alphaP, const double &z);
	Vector calcSourceAmplIn(const double &z, const int &pol);
	Vector calcIncAmplOut(const int &pol);
	Vector calcIncAmplIn(const double &z, const int &pol);
	Vector setSourceIntAmpl(const int &pol, const Complex &aP,
			const Complex &aM, const Complex &betaP, const double &alphaP,
			const double &zs, const Vector3d &pos);
	Vector3d setSourceIntFields(const int &pol, const Complex &aP,
			const Complex &aM, const Complex &betaP, const double &alphaP,
			const double &zs, const Vector3d &pos);
	Vector3d sourceExtFields(const Complex &e0, const Complex &beta,
			const double &alpha, const int &pol, const Vector3d &pos);
	Vector sourceExtAmpl(const Complex &e0, const Complex &beta,
			const double &alpha, const int &pol, const Vector3d &pos);
	void setLambda(const double &lambd);
	Complex getAngle(const int &l1, const Complex &beta1, const double &alpha1,
			const int &l2);
	int getMedium(const double &z);
	Matrix calcSMatrix(const int &pol);
	Matrix calcSMatrixU(const int &pol, const double &z);
	Matrix calcSMatrixL(const int &pol, const double &z);
	Matrix calcSMatrixInt(const int &pol, const double &z1, const double &z2);
	Matrix calcSl(const int &l, const int &u, const int &pol);
	Matrix calcSh(const int &i, const double &z);
	Vector calcFresnelCoefsTM(const int &i, const int &j);
	Vector calcFresnelCoefsTE(const int &i, const int &j);
	Complex calcKz(const int &lay, const bool &dirP);
	int nbrMediums;
	int *mediums;
	RVector Zi;
	Vector ni, ki, a0, aSa, aSb, aSPMa, aSPMb, a, aS;
	double lambda, alpha0;
	Complex beta0;
	Complex kx, ky, kz, Kzi, Kzj;
	//Matrix Sb,Sa,Sz,S;
	//Complex Ebout,Eaout;
	RefractiveIndex *ri;
};

#endif // MULTILAYER_H_INCLUDED
