#include "SpecialFunctions.h"
#include <iostream>
using namespace std;

double gammaln(const double &xx) {
	int j;
	double x, tmp, y, ser;
	double cof[14];
	cof[0] = 57.1562356658629235;
	cof[1] = -59.5979603554754912;
	cof[2] = 14.1360979747417471;
	cof[3] = -0.491913816097620199;
	cof[4] = .339946499848118887e-4;
	cof[5] = .465236289270485756e-4;
	cof[6] = -.983744753048795646e-4;
	cof[7] = .158088703224912494e-3;
	cof[8] = -.210264441724104883e-3;
	cof[9] = .217439618115212643e-3;
	cof[10] = -.164318106536763890e-3;
	cof[11] = .844182239838527433e-4;
	cof[12] = -.261908384015814087e-4;
	cof[13] = .368991826595316234e-5;
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j < 14; j++)
		ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005 * ser / x);
}

Complex fresnelCoefRP(const Complex &m, const Complex &cosBeta1) {
	Complex cosBeta2, sinBeta1, sinBeta2;
	sinBeta1 = sqrt(1. - cosBeta1 * cosBeta1);
	sinBeta2 = sinBeta1 / m;
	cosBeta2 = sqrt(1. - sinBeta2 * sinBeta2);
	if (imag(cosBeta2 * m) < 0.)
		cosBeta2 = -cosBeta2;
	if ((m * cosBeta1 - cosBeta2) < 1e-14)
		return 0.;
	return (m * cosBeta1 - cosBeta2) / (m * cosBeta1 + cosBeta2);
}

Complex fresnelCoefRS(const Complex &m, const Complex &cosBeta1) {
	Complex sinBeta1, cosBeta2, sinBeta2;
	sinBeta1 = sqrt(1. - cosBeta1 * cosBeta1);
	sinBeta2 = sinBeta1 / m;
	cosBeta2 = sqrt(1. - sinBeta2 * sinBeta2);
	if (imag(cosBeta2 * m) < 0.)
		cosBeta2 = -cosBeta2;
	if ((cosBeta1 - m * cosBeta2) <= 1e-14)
		return 0.;
	return (cosBeta1 - m * cosBeta2) / (cosBeta1 + m * cosBeta2);
}

Complex fresnelCoefTP(const Complex &m, const Complex &cosBeta1) {
	Complex cosBeta2, sinBeta1, sinBeta2;
	sinBeta1 = sqrt(1. - cosBeta1 * cosBeta1);
	sinBeta2 = sinBeta1 / m;
	cosBeta2 = sqrt(1. - sinBeta2 * sinBeta2);
	if (imag(cosBeta2) < 0.)
		cosBeta2 = -cosBeta2;
	if ((cosBeta1 + m * cosBeta2) < 1e-14)
		return 0.;
	return (2. * m * cosBeta1) / (cosBeta1 + m * cosBeta2);
}

Complex fresnelCoefTS(const Complex &m, const Complex &cosBeta1) {
	Complex sinBeta1, cosBeta2, sinBeta2;
	sinBeta1 = sqrt(1. - cosBeta1 * cosBeta1);
	sinBeta2 = sinBeta1 / m;
	cosBeta2 = sqrt(1. - sinBeta2 * sinBeta2);
	if (imag(cosBeta2) < 0.)
		cosBeta2 = -cosBeta2;
	if ((m * cosBeta1 + cosBeta2) <= 1e-14)
		return 0.;
	return (2. * m * cosBeta1) / (m * cosBeta1 + cosBeta2);
}

double gegenbauer(const double &alph, const int &n, const double &x) {
	if (n == 0)
		return 1.;
	if (n == 1)
		return 2. * alph * x;
	double Cn = 0., Cnm1, Cnm2;
	Cnm2 = 1.;
	Cnm1 = 2. * alph * x;
	for (int i = 2; i <= n; i++) {
		Cn = (2. * x * (double(i - 1) + alph) * Cnm1
				- (double(i - 2) + 2. * alph) * Cnm2) / double(i);
		Cnm2 = Cnm1;
		Cnm1 = Cn;
	}
	return Cn;
}

Complex Amnplus(const int &mm, const int &n, const double &vc) {
	int m = abs(mm);
	if (m > n)
		return 0.;
	Complex Ap = powi(n + m)
			* sqrt(double(n + n + 1) * fact(n - m) / (piDbl * fact(n + m)))
			* dblfact(2 * m - 1) / pow(vc / sqrt(1. - vc * vc), m)
			* gegenbauer(m + 0.5, n - m, 1. / vc) / vc;
	if (mm < 0)
		return powm1(m) * Ap;
	return Ap;
}

Complex Bmnplus(const int &mm, const int &n, const double &vc) {
	int m = abs(mm);
	Complex Bp = Amnplus(m + 1, n, vc) * sqrt(double(n + m + 1) * double(n - m))
			- Amnplus(m - 1, n, vc) * sqrt(double(n - m + 1) * double(n + m));
	if (mm < 0)
		return powm1(m) * Bp;
	return Bp;
}
