#ifndef VSWF_H_INCLUDED
#define VSWF_H_INCLUDED

#include "BesselFunctions.h"
#include "LegendreFunctions.h"
#include "Vector3d.h"
#include "Multilayer.h"

class VSWF {
public:
	VSWF();
	VSWF(const int &NN, const Complex &KR, const double &ct, const double &p);
	void setParameters(const int &N, const Complex &KR, const double &ct,
			const double &p);
	void calcVSWFs(const int &m, const int &n, const int &i,
			const bool &calcExp);
	int N;
	Complex kr, mult;
	double theta, phi, cost;
	SphericalBessel bes;
	LegendreFunctions leg;
	Vector3d Mmn, Nmn;
};

class VSWFintegral {
public:
	VSWFintegral();
	VSWFintegral(const int &NN, const Complex &KR, const double &ct,
			const double &p);
	void setParameters(const int &N, const Complex &KR, const double &ct,
			const double &p);
	void calcVSWFs(const int &m, const int &n, const int &i,
			const bool &calcExp);
	int N, Ng, Nga, Ngb;
	Complex kr, mult;
	double theta, phi, cost;
	SphericalBessel bes;
	LegendreFunctionsC *legA, *legB;
	Vector3d Mmn, Nmn;
	RMatrix gaussLegA, gaussLegB;
};

class VSWFsub {
public:
	VSWFsub();
	VSWFsub(const int &NN, const Complex &KM1, const Complex &KM2,
			const double &X, const double &Xj, const double &Y,
			const double &Yj, const double &Z, const double &Zj,
			const double &zsub);
	void setParameters(const int &NN, const Complex &KM1, const Complex &KM2,
			const double &X, const double &Xj, const double &Y,
			const double &Yj, const double &Z, const double &Zj,
			const double &zsub);
	void calcVSWFsRefl(const int &m, const int &n);
	void calcVSWFsTrans(const int &m, const int &n);
	void calcVSWFsIntegral(const int &m, const int &n);
	int N, Ng, Nga, Ngb;
	Complex km1, km2, mult;
	double theta, phi, r;
	double x, xj, y, yj, z, zj, zSub;
	bool reflection;
	Bessel bes;
	LegendreFunctionsC *leg;
	LegendreFunctionsC *legR, *legT;
	Vector3d Mmn, Nmn;
	RMatrix gaussLag, gaussLeg;
};

class VSWFsubinf {
public:
	VSWFsubinf();
	VSWFsubinf(const int &NN, const Complex &KM1, const Complex &KM2,
			const double &t, const double &p, const double &Zj);
	void setParameters(const int &NN, const Complex &KM1, const Complex &KM2,
			const double &t, const double &p, const double &Zj);
	void calcVSWFsRefl(const int &m, const int &n);
	void calcVSWFsTrans(const int &m, const int &n);
	int N;
	Complex km1, km2, mult, rs, rp, ts, tp;
	double theta, phi, r, zj, cost;
	Complex cosT, sinT;
	Vector3d Mmn, Nmn;
	LegendreFunctions legR, legT;
};

class VSWFinf {
public:
	VSWFinf();
	VSWFinf(const int &NN, const Complex &t, const double &p);
	void setParameters(const int &N, const Complex &t, const double &p);
	void calcVSWFs(const int &m, const int &n);
	int N;
	Complex mult, theta, cost;
	double phi;
	LegendreFunctionsC leg;
	Vector3d Mmn, Nmn;
};

class VSWFds {
public:
	VSWFds();
	VSWFds(const int &NN, const Complex &K, const double &rr,
			const Vector &Sources, const double &ct, const double &p);
	void setParameters(const int &NN, const Complex &K, const double &rr,
			const Vector &Sources, const double &ct, const double &p);
	void calcVSWFsds(const int &m, const int &n, const int &i,
			const bool &calcExp);
	Complex k, rho;
	double theta, phi, cost, sint, z, r;
	int N, m;
	Vector3d Mmn, Nmn;
	Vector kr;
	Vector sources, R, costC, sintC, thetaC;
	SphericalBessel *bes;
	LegendreFunctionsC *leg;
	VSWF vswf;
};

Vector3d calcMmnds(const int &m, const int &i, const Complex &k,
		const double &r, const Complex &zn, const double &ct);
Vector3d calcNmnds(const int &m, const int &i, const Complex &k,
		const double &r, const Complex &zn, const double &ct);

#endif // VSWF_H_INCLUDED
