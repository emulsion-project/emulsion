#ifndef MULTILAYER1_H_INCLUDED
#define MULTILAYER1_H_INCLUDED

#include "Matrix.h"
#include "RealContainers.h"
#include <fstream>
#include "RefractiveIndex.h"
#include "Vector3d.h"

class Multilayer1 {
public:
	Multilayer1();
	Multilayer1(const string &fname, const double &lambd);
	void setParameters(const string &fname, const double &lambd);
	void setLambda(const double &lambd);
	void setKxy(const Complex &k);
	int getMedium(const double &z);

	Vector sourceExtOut(const Vector &a, const int &pol);
	Vector sourceExtIn(const Vector &a, const int &pol, const double &z);
	Vector sourceInt(const int &pol, const Vector &a, const double &zs,
			const double &z);

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
	double lambda;
	Complex kxy, kz, Kzi, Kzj;
	//Matrix Sb,Sa,Sz,S;
	//Complex Ebout,Eaout;
	RefractiveIndex *ri;
};

#endif // MULTILAYER1_H_INCLUDED
