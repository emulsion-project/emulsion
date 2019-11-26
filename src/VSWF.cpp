#include "VSWF.h"
#include <iostream>
using namespace std;

VSWF::VSWF() {
}

VSWF::VSWF(const int &NN, const Complex &KR, const double &ct,
		const double &p) {
	N = NN;
	kr = KR;
	cost = ct;
	phi = p;
	bes.assign(N, kr);
	leg.assign(N, cost);
}

void VSWF::setParameters(const int &NN, const Complex &KR, const double &ct,
		const double &p) {
	N = NN;
	kr = KR;
	cost = ct;
	phi = p;
	bes.assign(N, kr);
	leg.assign(N, cost);
}

void VSWF::calcVSWFs(const int &m, const int &n, const int &i,
		const bool &calcExp) {
	mult = 1. / sqrt(2. * double(n * n + n));
	if (calcExp)
		mult *= exp(i_ * double(m) * phi);
	Mmn.t[0] = 0.;
	if (i == 1) {
		Mmn.t[1] = mult * bes.jn(n) * i_ * double(m) * leg.Pimn(n, m);
		Mmn.t[2] = -mult * bes.jn(n) * leg.Taumn(n, m);
		Nmn.t[0] = mult * bes.jn(n) / kr * double(n * n + n) * leg.Pmn(n, m);
		Nmn.t[1] = mult * bes.djn(n) / kr * leg.Taumn(n, m);
		Nmn.t[2] = mult * bes.djn(n) / kr * i_ * double(m) * leg.Pimn(n, m);
	} else if (i == 3) {
		Mmn.t[1] = mult * bes.h1n(n) * i_ * double(m) * leg.Pimn(n, m);
		Mmn.t[2] = -mult * bes.h1n(n) * leg.Taumn(n, m);
		Nmn.t[0] = mult * bes.h1n(n) / kr * double(n * n + n) * leg.Pmn(n, m);
		Nmn.t[1] = mult * bes.dh1n(n) / kr * leg.Taumn(n, m);
		Nmn.t[2] = mult * bes.dh1n(n) / kr * i_ * double(m) * leg.Pimn(n, m);
	}
	return;
}

/*VSWFintegral::VSWFintegral() {}

 VSWFintegral::VSWFintegral(const int &NN, const Complex &KR, const double &ct, const double &p) {
 N = NN;
 kr = KR;
 cost = ct;
 phi = p;

 Ng = 25;
 Nga = Ngb = 25;
 legA = new LegendreFunctionsC[Nga];
 legB = new LegendreFunctionsC[Ngb];
 for (int i=0 ; i<Nga ; i++) {
 legA[i].assign(N,t);
 }
 for (int i=0 ; i<Ngb ; i++) {
 legB[i].assign(N,t);
 }
 }

 void VSWFintegral::setParameters(const int &NN, const Complex &KR, const double &ct, const double &p) {
 N = NN;
 kr = KR;
 cost = ct;
 phi = p;

 Ng = 25;
 Nga = Ngb = 25;
 legA = new LegendreFunctionsC[Nga];
 legB = new LegendreFunctionsC[Ngb];
 for (int i=0 ; i<Nga ; i++) {
 legA[i].assign(N,t);
 }
 for (int i=0 ; i<Ngb ; i++) {
 legB[i].assign(N,t);
 }
 }

 void VSWFintegral::calcVSWFs(const int &m, const int &n, const int &i, const bool &calcExp) {
 Mmn = Nmn = Vector3d(0.);
 Vector3d MM(0.),NN(0.);
 Complex fact,t,mult,tp,ts,q;
 double alpha;
 for (int ka=0 ; ka<Ng ; ka++) {
 alpha = gaussLegA.t[ka][0];
 for (int kb=0 ; kb<Ngb ; kb++) {
 t = gaussLegB.t[kb][0];
 mult = exp(i_*double(m)*alpha);
 mult *= exp(i_*km1*r*(cos(alpha)+sin(alpha))*sqrt(1.-t*t));
 mult *= gaussLegA.t[ka][1]*gaussLegB.t[kb][1];
 mult *= -1./(2.*piDbl*powi(n+1)*sqrt(2*n*(n+1)));
 Mmn.t[0] += mult*(double(m)*leg.Pimn(n,m)*t*cos(alpha)-i_*leg.Taumn(n,m)*sin(alpha));
 Mmn.t[1] += mult*(double(m)*leg.Pimn(n,m)*t*sin(alpha)+i_*leg.Taumn(n,m)*cos(alpha));
 Mmn.t[2] -= mult*double(m)*leg.Pimn(n,m)*sqrt(1.-t*t);

 Nmn.t[0] += mult*(leg.Taumn(n,m)*t*cos(alpha)-i_*double(m)*leg.Pimn(n,m)*sin(alpha));
 Nmn.t[1] += mult*(leg.Taumn(n,m)*t*sin(alpha)+i_*double(m)*leg.Pimn(n,m)*cos(alpha));
 Nmn.t[2] -= mult*leg.Taumn(n,m)*sqrt(1.-t*t);
 }
 }
 return;
 }*/

VSWFsub::VSWFsub() {
}

VSWFsub::VSWFsub(const int &NN, const Complex &KM1, const Complex &KM2,
		const double &X, const double &Xj, const double &Y, const double &Yj,
		const double &Z, const double &Zj, const double &zsub) {
	N = NN;
	km1 = KM1;
	km2 = KM2;
	x = X;
	xj = Xj;
	y = Y;
	yj = Yj;
	z = Z;
	zj = Zj, zSub = zsub;
	r = sqrt((x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj));
	theta = acos(z - zj) / r;
	phi = atan2(x - xj, y - yj);
	Ng = 25;
	Nga = Ngb = 25;
	gaussLag = gaussLagQuad(Ngb, 0.);
	gaussLeg = gaussLegQuad(Nga, 0., 2. * piDbl);
	legR = new LegendreFunctionsC[Ngb];
	legT = new LegendreFunctionsC[Ngb];
	Complex q, t;
	q = i_ * km1 * (2. * (zSub - zj) + (zj - z));
	for (int i = 0; i < Ngb; i++) {
		t = 1. - gaussLag.t[i][0] / q;
		legR[i].assign(N, t);
	}
	q = i_ * km1 * (zSub - zj);
	for (int i = 0; i < Ngb; i++) {
		t = 1. - gaussLag.t[i][0] / q;
		legT[i].assign(N, t);
	}
}

void VSWFsub::setParameters(const int &NN, const Complex &KM1,
		const Complex &KM2, const double &X, const double &Xj, const double &Y,
		const double &Yj, const double &Z, const double &Zj,
		const double &zsub) {
	N = NN;
	km1 = KM1;
	km2 = KM2;
	x = X;
	xj = Xj;
	y = Y;
	yj = Yj;
	z = Z;
	zj = Zj, zSub = zsub;
	r = sqrt((x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj));
	theta = acos(z - zj) / r;
	phi = atan2(x - xj, y - yj);
	Ng = 25;
	Nga = Ngb = 25;
	gaussLag = gaussLagQuad(Ngb, 0.);
	gaussLeg = gaussLegQuad(Nga, 0., 2. * piDbl);
	leg = new LegendreFunctionsC[Ngb];
	legR = new LegendreFunctionsC[Ngb];
	legT = new LegendreFunctionsC[Ngb];
	Complex q, t;
	q = i_ * km1 * (2. * (zSub - zj) + (zj - z));
	for (int i = 0; i < Ngb; i++) {
		t = 1. - gaussLag.t[i][0] / q;
		legR[i].assign(N, t);
	}
	q = i_ * km1 * (zSub - zj);
	for (int i = 0; i < Ngb; i++) {
		t = 1. - gaussLag.t[i][0] / q;
		legT[i].assign(N, t);
	}
}

void VSWFsub::calcVSWFsRefl(const int &m, const int &n) {
	Mmn = Nmn = Vector3d(0.);
	Vector3d MM(0.), NN(0.);
	Complex fact, t, mult, rp, rs, q;
	double alpha;
	q = i_ * km1 * (2. * (zSub - zj) + (zj - z));
	for (int kb = 0; kb < Ngb; kb++) {
		t = 1. - gaussLag.t[kb][0] / q;
		rp = fresnelCoefRP(km2 / km1, t);
		rs = fresnelCoefRS(km2 / km1, t);
		for (int ka = 0; ka < Nga; ka++) {
			alpha = gaussLeg.t[ka][0];
			mult = exp(i_ * double(m) * alpha);
			mult *= exp(
					i_ * km1 * ((x - xj) * cos(alpha) + (y - yj) * sin(alpha))
							* sqrt(1. - t * t));
			mult *= gaussLeg.t[ka][1] * gaussLag.t[kb][1] * exp(q) / q;
			mult *= -1. / (2. * piDbl * powi(n + 1) * sqrt(2 * n * (n + 1)));
			Mmn.t[0] += mult
					* (-double(m) * legR[kb].Pimn(n, m) * t * cos(alpha) * rp
							- i_ * legR[kb].Taumn(n, m) * sin(alpha) * rs);
			Mmn.t[1] += mult
					* (double(m) * legR[kb].Pimn(n, m) * t * sin(alpha) * rp
							+ i_ * legR[kb].Taumn(n, m) * cos(alpha) * rs);
			Mmn.t[2] -= mult * double(m) * legR[kb].Pimn(n, m)
					* sqrt(1. - t * t) * rp;

			Nmn.t[0] += mult
					* (legR[kb].Taumn(n, m) * t * cos(alpha) * rp
							- i_ * double(m) * legR[kb].Pimn(n, m) * sin(alpha)
									* rs);
			Nmn.t[1] += mult
					* (legR[kb].Taumn(n, m) * t * sin(alpha) * rp
							+ i_ * double(m) * legR[kb].Pimn(n, m) * cos(alpha)
									* rs);
			Nmn.t[2] -= mult * legR[kb].Taumn(n, m) * sqrt(1. - t * t) * rp;
		}
	}
	return;
}

void VSWFsub::calcVSWFsTrans(const int &m, const int &n) {
	Mmn = Nmn = Vector3d(0.);
	Vector3d MM(0.), NN(0.);
	Complex fact, t, mult, tp, ts, q, cosT, sinT;
	double alpha;
	q = i_ * km1 * (zSub - zj);
	for (int kb = 0; kb < Ngb; kb++) {
		t = 1. - gaussLag.t[kb][0] / q;
		sinT = sqrt(1. - t * t) * km1 / km2;
		cosT = sqrt(1. - sinT * sinT);
		if (cosT.im < 0.)
			cosT *= -1.;
		tp = fresnelCoefTP(km2 / km1, t);
		ts = fresnelCoefTS(km2 / km1, t);
		for (int ka = 0; ka < Nga; ka++) {
			alpha = gaussLeg.t[ka][0];
			mult = exp(i_ * double(m) * alpha);
			mult *= exp(
					i_ * km1 * ((x - xj) * cos(alpha) + (y - yj) * sin(alpha))
							* sqrt(1. - t * t));
			mult *= exp(i_ * km2 * (z - zSub) * cosT);
			mult *= gaussLeg.t[ka][1] * gaussLag.t[kb][1] * exp(q) / q;
			mult *=
					-1.
							/ (2. * piDbl * powi(n + 1)
									* sqrt(double(2 * n * (n + 1))));
			Mmn.t[0] += mult
					* (-double(m) * legT[kb].Pimn(n, m) * t * cos(alpha) * tp
							- i_ * legT[kb].Taumn(n, m) * sin(alpha) * ts);
			Mmn.t[1] += mult
					* (double(m) * legT[kb].Pimn(n, m) * t * sin(alpha) * tp
							+ i_ * legT[kb].Taumn(n, m) * cos(alpha) * ts);
			Mmn.t[2] -= mult * double(m) * legT[kb].Pimn(n, m)
					* sqrt(1. - t * t) * tp;

			Nmn.t[0] += mult
					* (legT[kb].Taumn(n, m) * t * cos(alpha) * tp
							- i_ * double(m) * legT[kb].Pimn(n, m) * sin(alpha)
									* ts);
			Nmn.t[1] += mult
					* (legT[kb].Taumn(n, m) * t * sin(alpha) * tp
							+ i_ * double(m) * legT[kb].Pimn(n, m) * cos(alpha)
									* ts);
			Nmn.t[2] -= mult * legT[kb].Taumn(n, m) * sqrt(1. - t * t) * tp;
		}
	}
	return;
}

void VSWFsub::calcVSWFsIntegral(const int &m, const int &n) {
	/*Mmn = Nmn = Vector3d(0.);
	 Vector3d MM(0.),NN(0.);
	 Complex fact,t,mult,tp,ts,q;
	 double alpha;
	 q = i_*km1*(z-zj);
	 for (int ka=0 ; ka<Nga ; ka++) {
	 alpha = gaussLeg.t[ka][0];
	 for (int kb=0 ; kb<Ngb ; kb++) {
	 t = 1.-gaussLag.t[kb][0]/q;
	 //leg.assign(n,t);
	 mult = exp(i_*double(m)*alpha);
	 mult *= exp(i_*km1*((x-xj)*cos(alpha)+(y-yj)*sin(alpha))*sqrt(1.-t*t));
	 mult *= gaussLeg.t[ka][1]*gaussLag.t[kb][1]*exp(q)/q;
	 mult *= -1./(2.*piDbl*powi(n+1)*sqrt(2*n*(n+1)));
	 Mmn.t[0] += mult*(double(m)*leg.Pimn(n,m)*t*cos(alpha)-i_*leg.Taumn(n,m)*sin(alpha));
	 Mmn.t[1] += mult*(double(m)*leg.Pimn(n,m)*t*sin(alpha)+i_*leg.Taumn(n,m)*cos(alpha));
	 Mmn.t[2] -= mult*double(m)*leg.Pimn(n,m)*sqrt(1.-t*t);

	 Nmn.t[0] += mult*(leg.Taumn(n,m)*t*cos(alpha)-i_*double(m)*leg.Pimn(n,m)*sin(alpha));
	 Nmn.t[1] += mult*(leg.Taumn(n,m)*t*sin(alpha)+i_*double(m)*leg.Pimn(n,m)*cos(alpha));
	 Nmn.t[2] -= mult*leg.Taumn(n,m)*sqrt(1.-t*t);
	 }
	 }*/
	return;
}

VSWFsubinf::VSWFsubinf() {
}

VSWFsubinf::VSWFsubinf(const int &NN, const Complex &KM1, const Complex &KM2,
		const double &t, const double &p, const double &Zj) {
	N = NN;
	km1 = KM1;
	km2 = KM2;
	theta = t;
	cost = cos(theta);
	sinT = sqrt(1. - cost * cost) * km2 / km1;
	cosT = sqrt(1. - sinT * sinT);
	phi = p;
	zj = Zj;
	sinT = sqrt(1. - cost * cost) * km1 / km2;
	cosT = sqrt(1. - sinT * sinT);
	legR.assign(N, -cost);
	legT.assign(N, cost);
}

void VSWFsubinf::setParameters(const int &NN, const Complex &KM1,
		const Complex &KM2, const double &t, const double &p,
		const double &Zj) {
	N = NN;
	km1 = KM1;
	km2 = KM2;
	theta = t;
	phi = p;
	zj = Zj;
	sinT = sqrt(1. - cost * cost) * km1 / km2;
	cosT = sqrt(1. - sinT * sinT);
	if (cosT.im < 0.)
		cosT *= -1.;
	legR.assign(N, -cost);
	legR.assign(N, cost);
}

void VSWFsubinf::calcVSWFsRefl(const int &m, const int &n) {
	mult = 1. / sqrt(2. * double(n * n + n));
	mult *= exp(i_ * double(m) * phi);
	mult *= exp(-2. * i_ * km1 * zj * cost);
	rp = fresnelCoefRP(km2 / km1, -cost);
	rs = fresnelCoefRS(km2 / km1, -cost);
	Mmn.t[0] = Nmn.t[0] = 0.;
	Mmn.t[1] = mult * i_ * double(m) * legR.Pimn(n, m) * rp;
	Mmn.t[2] = -mult * legR.Taumn(n, m) * rs;
	Nmn.t[1] = mult * legR.Taumn(n, m) * rp;
	Nmn.t[2] = mult * i_ * double(m) * legR.Pimn(n, m) * rs;
	return;
}

void VSWFsubinf::calcVSWFsTrans(const int &m, const int &n) {

	mult = 1. / sqrt(2. * double(n * n + n));
	mult *= exp(i_ * double(m) * phi);
	mult *= exp(i_ * zj * (km2 * cost - km1 * cosT));
	tp = fresnelCoefTP(km2 / km1, cosT);
	ts = fresnelCoefTS(km2 / km1, cosT);
	Mmn.t[0] = Nmn.t[0] = 0.;
	Mmn.t[1] = mult * i_ * double(m) * legT.Pimn(n, m) * tp;
	Mmn.t[2] = -mult * legT.Taumn(n, m) * ts;
	Nmn.t[1] = mult * legT.Taumn(n, m) * tp;
	Nmn.t[2] = mult * i_ * double(m) * legT.Pimn(n, m) * ts;
	return;
}

VSWFinf::VSWFinf() {
}

VSWFinf::VSWFinf(const int &NN, const Complex &t, const double &p) {
	N = NN;
	theta = t;
	cost = cos(theta);
	phi = p;
	leg.assign(N, cost);
}

void VSWFinf::setParameters(const int &NN, const Complex &t, const double &p) {
	N = NN;
	theta = t;
	cost = cos(theta);
	phi = p;
	leg.assign(N, cost);
}

void VSWFinf::calcVSWFs(const int &m, const int &n) {
	mult = 1. / sqrt(2. * double(n * n + n));
	mult *= exp(i_ * double(m) * phi);
	Mmn.t[0] = Nmn.t[0] = 0.;
	Mmn.t[1] = mult * i_ * double(m) * leg.Pimn(n, m);
	Mmn.t[2] = -mult * leg.Taumn(n, m);
	Nmn.t[1] = mult * leg.Taumn(n, m);
	Nmn.t[2] = mult * i_ * double(m) * leg.Pimn(n, m);
	return;
}

VSWFds::VSWFds() {
	bes = NULL;
	leg = NULL;
}

VSWFds::VSWFds(const int &NN, const Complex &K, const double &rr,
		const Vector &Sources, const double &ct, const double &p) {
	N = NN;
	k = K;
	r = rr;
	sources = Sources;
	theta = ct;
	cost = cos(theta);
	sint = sin(theta);
	z = r * cost;
	rho = r * sint;
	phi = p;
	bes = new SphericalBessel[N];
	leg = new LegendreFunctionsC[N];
	kr = R = costC = sintC = Vector(N, 0.);
	for (int i = 0; i < N; i++) {
		R.t[i] = sqrt(rho * rho + (z - sources.t[i]) * (z - sources.t[i]));
		//if ((R.t[i]) < 1e-14) R.t[i] = 1e-14;
		kr.t[i] = k * R.t[i];
		costC.t[i] = (z - sources.t[i]) / R.t[i];
		sintC.t[i] = rho / R.t[i];
		bes[i].assign(N, kr.t[i]);
		leg[i].assign(N, costC.t[i]);
	}
}

void VSWFds::setParameters(const int &NN, const Complex &K, const double &rr,
		const Vector &Sources, const double &ct, const double &p) {
	N = NN;
	k = K;
	r = rr;
	sources = Sources;
	theta = ct;
	cost = cos(theta);
	sint = sin(theta);
	//sint = sqrt(1.-cost*cost);
	z = r * cost;
	rho = r * sint;
	phi = p;
	if (bes != NULL)
		delete[] bes;
	if (leg != NULL)
		delete[] leg;
	bes = new SphericalBessel[N];
	leg = new LegendreFunctionsC[N];
	kr = R = costC = sintC = Vector(N, 0.);
	for (int i = 0; i < N; i++) {
		R.t[i] = sqrt(rho * rho + (z - sources.t[i]) * (z - sources.t[i]));
		//if ((R.t[i]) < 1e-14) R.t[i] = 1e-14;
		kr.t[i] = k * R.t[i];
		costC.t[i] = (z - sources.t[i]) / R.t[i];
		sintC.t[i] = rho / R.t[i];
		bes[i].assign(N, kr.t[i]);
		leg[i].assign(N, costC.t[i]);
	}
}

void VSWFds::calcVSWFsds(const int &m, const int &n, const int &i,
		const bool &calcExp) {
	int nn = abs(m);
	if (nn == 0)
		nn = 1;
	Complex mult = 1. / sqrt(2. * double(nn * nn + nn));
	if (calcExp)
		mult *= exp(i_ * double(m) * phi);
	Complex stC = sint * costC.t[n - 1] - cost * sintC.t[n - 1];
	Complex ctC = cost * costC.t[n - 1] + sint * sintC.t[n - 1];
	Complex p = double(nn * nn + nn) * leg[n - 1].Pmn(nn, m);
	Complex pi = i_ * double(m) * leg[n - 1].Pimn(nn, m);
	Complex tau = leg[n - 1].Taumn(nn, m);
	Complex jh, djh;
	if (i == 1) {
		jh = bes[n - 1].jn(nn);
		djh = bes[n - 1].djn(nn);
	} else if (i == 3) {
		jh = bes[n - 1].h1n(nn);
		djh = bes[n - 1].dh1n(nn);
	}
	Mmn.t[0] = mult * jh * pi * stC;
	Mmn.t[1] = mult * jh * pi * ctC;
	Mmn.t[2] = -mult * jh * tau;
	Nmn.t[0] = mult * (jh * p * ctC + djh * tau * stC) / kr.t[n - 1];
	Nmn.t[1] = mult * (-jh * p * stC + djh * tau * ctC) / kr.t[n - 1];
	Nmn.t[2] = mult * djh * pi / kr.t[n - 1];
	return;
}

Vector3d calcMmnds(const int &m, const int &i, const Complex &k,
		const double &r, const Complex &zn, const double &ct) {
	int n = abs(m);
	if (m == 0)
		n = 1;
	int NBes = n;
	double mr = double(m);
	double nr = double(n * n + n);
	double nm = 1. / sqrt(2. * nr);
	double cth = cos(ct);
	double sth = sin(ct);
	double ro = r * sth;
	double z = r * cth;
	Complex dz = z - zn;
	Complex RR = sqrt(ro * ro + dz * dz);
	Complex sint = ro / RR;
	Complex cost = dz / RR;
	Complex argJ = k * RR;
	SphericalBessel bes(NBes, argJ);
	LegendreFunctionsC leg(n, cost);
	Complex sinc = sth * cost - cth * sint;
	Complex cosc = cth * cost + sth * sint;
	Complex fp = i_ * mr * leg.Pimn(n, m) * nm;
	Complex ft = leg.Taumn(n, m) * nm;
	Complex jh;
	if (i == 1)
		jh = bes.jn(n);
	else
		jh = bes.h1n(n);
	Complex factp = jh * fp;
	Vector3d M(0.);
	M.t[0] = factp * sinc;
	M.t[1] = factp * cosc;
	M.t[2] = -jh * ft;
	return M;
}

Vector3d calcNmnds(const int &m, const int &i, const Complex &k,
		const double &r, const Complex &zn, const double &ct) {
	int n = abs(m);
	if (m == 0)
		n = 1;
	int NBes = n;
	double mr = double(m);
	double nr = double(n * n + n);
	double nm = 1. / sqrt(2. * nr);
	double cth = cos(ct);
	double sth = sin(ct);
	double ro = r * sth;
	double z = r * cth;
	Complex dz = z - zn;
	Complex RR = sqrt(ro * ro + dz * dz);
	Complex sint = ro / RR;
	Complex cost = dz / RR;
	Complex argJ = k * RR;
	SphericalBessel bes(NBes, argJ);
	LegendreFunctionsC leg(n, cost);
	Complex sinc = sth * cost - cth * sint;
	Complex cosc = cth * cost + sth * sint;
	Complex fp = i_ * mr * leg.Pimn(n, m) * nm;
	Complex ft = leg.Taumn(n, m) * nm;
	Complex fl = nr * leg.Pmn(n, m) * nm;
	Complex jh;
	if (i == 1)
		jh = bes.jn(n);
	else
		jh = bes.h1n(n);
	Complex djh;
	if (i == 1)
		djh = bes.djn(n);
	else
		djh = bes.dh1n(n);
	Complex factl = jh * fl;
	Complex factt = djh * ft;
	Vector3d M(0.);
	M.t[0] = (factl * cosc + factt * sinc) / argJ;
	M.t[1] = (-factl * sinc + factt * cosc) / argJ;
	M.t[2] = djh * fp / argJ;
	return M;
}
