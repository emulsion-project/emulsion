#include "IncidentRadiation.h"

#include <iostream>
using namespace std;

IncidentRadiation::IncidentRadiation() {
}

IncidentRadiation::IncidentRadiation(const int &inc, const double &lambd) {
	setParameters(inc, lambd);
}

void IncidentRadiation::setParameters(const int &inc, const double &lambd) {

	ml.setParameters("Multilayer/multilayerConfig.txt", lambd);

	N = 1;
	incidence = inc;
	lambda = lambd;
	ifstream ifs("Incidence/incidence.txt");
	int i = 1;
	string str;
	while (i < incidence) {
		getline(ifs, str);
		i++;
	}
	ifs >> type >> fileName;
	if (type == 0) { // plane wave
		ifstream ifs1(fileName.c_str());
		double re, im;
		ifs1 >> str >> beta0 >> str >> alpha0;
		beta0 *= piDbl / 180.;
		alpha0 *= piDbl / 180.;
		ifs1 >> str >> re >> im;
		Eb = Complex(re, im);
		ifs1 >> str >> re >> im;
		Ea = Complex(re, im);
		aincs = Vector(2, 0.);
		aincp = Vector(2, 0.);
		if (beta0 <= piDbl / 2.) {
			aincs.t[0] = Ea;
			aincp.t[0] = Eb;
			ml.setKxy(
					ml.ki.t[0]
							* sqrt(
									cos(alpha0) * sin(beta0) * cos(alpha0)
											* sin(beta0)
											+ sin(alpha0) * sin(beta0)
													* sin(alpha0)
													* sin(beta0)));
		} else {
			aincs.t[1] = Ea;
			aincp.t[1] = Eb;
			ml.setKxy(
					ml.ki.t[ml.nbrMediums - 1]
							* sqrt(
									cos(alpha0) * sin(beta0) * cos(alpha0)
											* sin(beta0)
											+ sin(alpha0) * sin(beta0)
													* sin(alpha0)
													* sin(beta0)));
		}
		aouts = ml.sourceExtOut(aincs, 0);
		aoutp = ml.sourceExtOut(aincp, 1);
	}
	if (type == 1) { // electron
		ifstream ifs1(fileName.c_str());
		ifs1 >> str >> E0 >> E >> str >> x0 >> y0 >> z0;
		E = 2. * piDbl * c * hbar / (lambd * e);
	}
}

void IncidentRadiation::setPlaneWaveCoeffs(const int &n, const Vector3d &pos) {
	N = n;
	int NM = (N + 1) * (N + 1) - 1;
	Amn = Vector(2 * NM);
	Complex mult, multn;

	int i;
	int lay = ml.getMedium(real(pos.t[2]));

	aincs = Vector(2, 0.);
	aincp = Vector(2, 0.);
	if (beta0 <= piDbl / 2.) {
		aincs.t[0] = Ea;
		aincp.t[0] = Eb;
		ml.setKxy(
				ml.ki.t[0]
						* sqrt(
								cos(alpha0) * sin(beta0) * cos(alpha0)
										* sin(beta0)
										+ sin(alpha0) * sin(beta0) * sin(alpha0)
												* sin(beta0)));
	} else {
		aincs.t[1] = Ea;
		aincp.t[1] = Eb;
		ml.setKxy(
				ml.ki.t[ml.nbrMediums - 1]
						* sqrt(
								cos(alpha0) * sin(beta0) * cos(alpha0)
										* sin(beta0)
										+ sin(alpha0) * sin(beta0) * sin(alpha0)
												* sin(beta0)));
	}
	aouts = ml.sourceExtOut(aincs, 0);
	aoutp = ml.sourceExtOut(aincp, 1);
	ap = ml.sourceExtIn(aincp, 1, pos.t[2].re);
	as = ml.sourceExtIn(aincs, 0, pos.t[2].re);
	ap = exp(i_ * ml.kxy * (pos.t[0] * cos(alpha0) + pos.t[1] * sin(alpha0)))
			* ap;
	as = exp(i_ * ml.kxy * (pos.t[0] * cos(alpha0) + pos.t[1] * sin(alpha0)))
			* as;
	LegendreFunctionsC leg(N, ml.calcKz(lay, true) / ml.ki.t[lay]);
	i = 0;

	for (int n = 1; n <= N; n++) {
		multn = -4. * powi(n) / sqrt(2. * double(n * n + n));
		for (int m = -n; m <= n; m++) {
			mult = multn * exp(-i_ * double(m) * alpha0);
			Amn.t[i] = mult
					* (i_ * double(m) * leg.Pimn(n, m) * ap.t[0]
							+ leg.Taumn(n, m) * as.t[0]);
			Amn.t[i + NM] = mult * i_
					* (leg.Taumn(n, m) * ap.t[0]
							- i_ * double(m) * leg.Pimn(n, m) * as.t[0]);
			i++;
		}
	}
	leg.assign(N, ml.calcKz(lay, false) / ml.ki.t[lay]);
	i = 0;
	for (int n = 1; n <= N; n++) {
		multn = -4. * powi(n) / sqrt(2. * double(n * n + n));
		for (int m = -n; m <= n; m++) {
			mult = multn * exp(-i_ * double(m) * alpha0);
			Amn.t[i] += mult
					* (i_ * double(m) * leg.Pimn(n, m) * ap.t[1]
							+ leg.Taumn(n, m) * as.t[1]);
			Amn.t[i + NM] += mult * i_
					* (leg.Taumn(n, m) * ap.t[1]
							- i_ * double(m) * leg.Pimn(n, m) * as.t[1]);
			i++;
		}
	}
	return;
}

void IncidentRadiation::setElectron(const int &n, const double &e0,
		const double &ee, const double &x, const double &y) {
	N = n;
	int NM = (N + 1) * (N + 1) - 1;

	E0 = e0;
	v0 = sqrt(2. * E0 * e / me);
	double vc = v0 / c;
	gamma = 1. / sqrt(1. - vc * vc);
	E = ee;
	omega = E * e / hbar;
	b = sqrt((x0 - x) * (x0 - x) + (y0 - y) * (y0 - y));
	phi0 = atan2(y0 - y, x0 - x);
	Amn = Vector(2 * NM);

	int i = 0;

	Complex mult;
	for (int n = 1; n <= N; n++) {
		for (int m = -n; m <= n; m++) {
			mult = -2. * piDbl * powi(1 - n) * omega
					/ (c * c * double(n * n + n))
					* besselKn(m, omega * b / (v0 * gamma))
					* exp(-i_ * double(m) * phi0 - i_ * omega * z0 / v0);
			Amn.t[i] = mult * 2. * vc * double(m) * Amnplus(m, n, vc);
			Amn.t[i + NM] = mult / gamma * Bmnplus(m, n, vc);
			i++;
		}
	}
}
