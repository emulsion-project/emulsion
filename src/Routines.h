#ifndef ROUTINES_H_INCLUDED
#define ROUTINES_H_INCLUDED

#include "MultipleScattering.h"
#include <iostream>
#include <fstream>
using namespace std;

void calcCrossSections(const string &fname, const double &lambda1,
		const double &lambda2, const double &step, const bool &normalize,
		const double &cstNorm);
void calcCrossSections(const string &fname, const double &lambda1,
		const double &lambda2, const double &step, const bool &normalize,
		const double &cstNorm, const int &nn, const int &pol);
void calcComplexCext(const string &fname, const double &lambda1,
		const double &lambda2, const double &step, const bool &normalize,
		const double &cstNorm);
void calcComplexCextN(const string &fname, const double &lambda1,
		const double &lambda2, const double &step, const bool &normalize,
		const double &cstNorm, const int &nn, const int &pol);
void calcCabsi(const string &fname, const double &lambda1,
		const double &lambda2, const double &step, const bool &normalize,
		const double &cstNorm);
void calcCextVSn(const string &fname, const int &N1, const int &N2,
		const int &pas, const bool &normalize, const double &cstNorm);
void calcEnergyLossProbability(const string &fname, const double &e1,
		const double &e2, const double &step);
void calcNearFieldPlane(const string &fname, const int &NT, const int &NP,
		const double &t1, const double &t2, const double &p1, const double &p2,
		const int &plane, const bool &cart);
void calcNearFieldParticleSurface(const string &fname, const int &NT,
		const int &NP, const bool &cart);
void calcFarFieldSurface(const string &fname, const int &N1, const int &N2);
void calcSurfaceCharges(const string &fname, const int &N1, const int &N2);
void calcEELSMapping(const string &fname, const int &NT, const int &NP,
		const double &t1, const double &t2, const double &p1, const double &p2);
Matrix polesAmpExtraction(const string &fname, const int &nPoles,
		const int &nPoints, const double &l1, const double &l2,
		const int &part);
Matrix polesAmpExtraction(MultipleScattering &ms, const int &nPoles,
		const int &nPoints, const double &l1, const double &l2,
		const int &part);
Vector poleAmpExtraction(const string &fname, const int &nPoints,
		const double &l1, const double &l2, const int &part);
Vector poleAmpExtraction(const string &fname, const double &z,
		const int &nPoints, const double &l1, const double &l2, const int &n,
		const int &m, const int &elec);

#endif // ROUTINES_H_INCLUDED
