#include "Transmission.h"
#include <iostream>
using namespace std;

Transmission::Transmission() {
	ml.setParameters("Multilayer/multilayerConfig.txt", 488e-9);
	N = 0;
	Ng = 0;
}

Transmission::Transmission(const double &lambd, const int &n, const int &ng,
		const double &xx, const double &yy, const double &zzi, const double zzj,
		const bool &self) {
	ml.setParameters("Multilayer/multilayerConfig.txt", lambd);
	N = 0;
	Ng = 0;
	setParameters(lambd, n, ng, xx, yy, zzi, zzj, self);
}

void Transmission::setParameters(const double &lambd, const int &n,
		const int &ng, const double &xx, const double &yy, const double &zzi,
		const double zzj, const bool &self) {
	ml.setLambda(lambd);
	lay = ml.getMedium(zi);
	xMax = 2.;
	if ((N != n) && (Ng != ng)) {
		Ng = ng;
		N = n;
		xInt = new Complex[Ng];
		leg = new LegendreFunctions[Ng];
		legC = new LegendreFunctionsC[Ng];
		gaussLeg = gaussLegQuad(Ng, 0., 1.);
		pasX = xMax / double(Ng - 1);
		for (int i = 0; i < Ng; i++) {
			xInt[i] = i_ * double(i) * pasX;
			leg[i].assign(N, gaussLeg.t[i][0]);
			legC[i].assign(N, i_ * xInt[i]);
		}
	}
	NM = (N + 1) * (N + 1) - 1;

	dx = xx;
	dy = yy;
	zi = zzi;
	zj = zzj;

	lay = ml.getMedium(zi);
	gaussLag = gaussLagQuad(Ng, 0.);
	if (self) {
		if (lay == ml.nbrMediums - 1)
			R = calcRSelf(0);
		else if (lay == 0)
			R = calcRSelf(1);
		else
			R = calcRSelf(0) + calcRSelf(1);
	} else {
		if (lay == ml.nbrMediums - 1)
			R = calcR(0);
		else if (lay == 0)
			R = calcR(1);
		else
			R = calcR(1);
	}
	return;
}

Matrix Transmission::calcRSelf(const bool &Up) {
	Matrix RR(2 * NM, 2 * NM, 0.);
	int line = 0, col = 0;
	Complex fact, t, mult, f, fpp, ftt, fpt, ftp, rp, rs, q;
	for (int k = 0; k < Ng; k++) {

		//ml.setSourceInt(1.,0.,1.,0.,acos(gaussLeg.t[k][0]),0.,zi);
		if (Up) {
			rp = ml.aSPMb.t[1];
			rs = ml.aSPMa.t[1];
		} else {
			rp = ml.aSPMb.t[0];
			rs = ml.aSPMa.t[0];
		}
		line = 0;
		if (Up)
			q = 2. * i_ * ml.ki.t[lay] * fabs(ml.Zi.t[lay] - zj)
					* gaussLeg.t[k][0];
		else
			q = 2. * i_ * ml.ki.t[lay - 1] * fabs(zj - ml.Zi.t[lay - 1])
					* gaussLeg.t[k][0];

		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				col = 0;
				for (int l = 1; l <= N; l++) {
					for (int kk = -l; kk <= l; kk++) {
						if (m == kk) {
							if (Up) {
								mult = 2. * gaussLeg.t[k][1] * powi(l - n)
										* powm1(m - l)
										/ sqrt(
												double(n * n + n)
														* double(l * l + l));
								if ((lay == 0) || (lay == ml.nbrMediums - 1))
									mult *= exp(q);

								RR.t[line][col] -= mult
										* (double(m * m) * leg[k].Pimn(n, m)
												* leg[k].Pimn(l, m) * rp
												- leg[k].Taumn(n, m)
														* leg[k].Taumn(l, m)
														* rs);
								RR.t[line][col + NM] -=
										mult * double(m)
												* (-leg[k].Pimn(n, m)
														* leg[k].Taumn(l, m)
														* rp
														+ leg[k].Taumn(n, m)
																* leg[k].Pimn(l,
																		m) * rs);
								RR.t[line + NM][col] -= mult * double(m)
										* (leg[k].Taumn(n, m)
												* leg[k].Pimn(l, m) * rp
												- leg[k].Pimn(n, m)
														* leg[k].Taumn(l, m)
														* rs);
								RR.t[line + NM][col + NM] -=
										mult
												* (-leg[k].Taumn(n, m)
														* leg[k].Taumn(l, m)
														* rp
														+ double(m * m)
																* leg[k].Pimn(n,
																		m)
																* leg[k].Pimn(l,
																		m) * rs);
							} else {
								mult = 2. * gaussLeg.t[k][1] * powi(l - n)
										* powm1(m - n)
										/ sqrt(
												double(n * n + n)
														* double(l * l + l));
								if ((lay == 0) || (lay == ml.nbrMediums - 1))
									mult *= exp(q);
								RR.t[line][col] -= mult
										* (double(m * m) * leg[k].Pimn(n, m)
												* leg[k].Pimn(l, m) * rp
												- leg[k].Taumn(n, m)
														* leg[k].Taumn(l, m)
														* rs);
								RR.t[line][col + NM] -=
										mult * double(m)
												* (leg[k].Pimn(n, m)
														* leg[k].Taumn(l, m)
														* rp
														- leg[k].Taumn(n, m)
																* leg[k].Pimn(l,
																		m) * rs);
								RR.t[line + NM][col] -= mult * double(m)
										* (-leg[k].Taumn(n, m)
												* leg[k].Pimn(l, m) * rp
												+ leg[k].Pimn(n, m)
														* leg[k].Taumn(l, m)
														* rs);
								RR.t[line + NM][col + NM] -=
										mult
												* (-leg[k].Taumn(n, m)
														* leg[k].Taumn(l, m)
														* rp
														+ double(m * m)
																* leg[k].Pimn(n,
																		m)
																* leg[k].Pimn(l,
																		m) * rs);
							}
						}
						col++;
					}
				}
				line++;
			}
		}
	}

	for (int k = 0; k < Ng; k++) {
//    ml.setSourceInt(1.,0.,1.,0.,acos(xInt[k]),0.,zi);
		if (Up) {
			rp = ml.aSPMb.t[1];
			rs = ml.aSPMa.t[1];
		} else {
			rp = ml.aSPMb.t[0];
			rs = ml.aSPMa.t[0];
		}
		//cout << xInt[k] << "\t" << Up << "\t" << rp << endl;
		line = 0;
		if (Up)
			q = 2. * i_ * ml.ki.t[lay] * fabs(ml.Zi.t[lay] - zj) * xInt[k];
		else
			q = 2. * i_ * ml.ki.t[lay - 1] * fabs(zj - ml.Zi.t[lay - 1])
					* xInt[k];
		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				col = 0;
				for (int l = 1; l <= N; l++) {
					for (int kk = -l; kk <= l; kk++) {
						if (m == kk) {
							if (Up) {
								mult = pasX * 2. * i_ * powi(l - n)
										* powm1(m - l)
										/ sqrt(
												double(n * n + n)
														* double(l * l + l));
								if ((k == 0) || (k == Ng - 1))
									mult *= 0.5;
								if ((lay == 0) || (lay == ml.nbrMediums - 1))
									mult *= exp(q);
								RR.t[line][col] += mult
										* (double(m * m) * legC[k].Pimn(n, m)
												* legC[k].Pimn(l, m) * rp
												- legC[k].Taumn(n, m)
														* legC[k].Taumn(l, m)
														* rs);
								RR.t[line][col + NM] += mult * double(m)
										* (-legC[k].Pimn(n, m)
												* legC[k].Taumn(l, m) * rp
												+ legC[k].Taumn(n, m)
														* legC[k].Pimn(l, m)
														* rs);
								RR.t[line + NM][col] += mult * double(m)
										* (legC[k].Taumn(n, m)
												* legC[k].Pimn(l, m) * rp
												- legC[k].Pimn(n, m)
														* legC[k].Taumn(l, m)
														* rs);
								RR.t[line + NM][col + NM] += mult
										* (-legC[k].Taumn(n, m)
												* legC[k].Taumn(l, m) * rp
												+ double(m * m)
														* legC[k].Pimn(n, m)
														* legC[k].Pimn(l, m)
														* rs);
								//if (n==1 && m==1 && l==1) cout << xInt[k] << "\t" << mult << "\t" << mult*(double(m*m)*legC[k].Pimn(n,m)*legC[k].Pimn(l,m)*rp-legC[k].Taumn(n,m)*legC[k].Taumn(l,m)*rs) << endl;
								//cout << legC[k].Pimn(n,m) << endl;
							} else {
								mult = pasX * 2. * powi(l - n) * i_
										* powm1(m - n)
										/ sqrt(
												double(n * n + n)
														* double(l * l + l));
								if ((k == 0) || (k == Ng - 1))
									mult *= 0.5;
								if ((lay == 0) || (lay == ml.nbrMediums - 1))
									mult *= exp(q);

								RR.t[line][col] += mult
										* (double(m * m) * legC[k].Pimn(n, m)
												* legC[k].Pimn(l, m) * rp
												- legC[k].Taumn(n, m)
														* legC[k].Taumn(l, m)
														* rs);
								RR.t[line][col + NM] += mult * double(m)
										* (legC[k].Pimn(n, m)
												* legC[k].Taumn(l, m) * rp
												- legC[k].Taumn(n, m)
														* legC[k].Pimn(l, m)
														* rs);
								RR.t[line + NM][col] += mult * double(m)
										* (-legC[k].Taumn(n, m)
												* legC[k].Pimn(l, m) * rp
												+ legC[k].Pimn(n, m)
														* legC[k].Taumn(l, m)
														* rs);
								RR.t[line + NM][col + NM] += mult
										* (-legC[k].Taumn(n, m)
												* legC[k].Taumn(l, m) * rp
												+ double(m * m)
														* legC[k].Pimn(n, m)
														* legC[k].Pimn(l, m)
														* rs);
							}
						}
						col++;
					}
				}
				line++;
			}
		}
	}

	return RR;
}

Matrix Transmission::calcR(const bool &Up) {
	Matrix RR(2 * NM, 2 * NM, 0.);
	int line = 0, col = 0;
	double rauij;
	Complex mult, rp, rs, q, eiphi;
	for (int k = 0; k < Ng; k++) {
//    ml.setSourceInt(1.,0.,1.,0.,acos(gaussLeg.t[k][0]),0.,zi);
		if (Up) {
			rp = ml.aSPMb.t[1];
			rs = ml.aSPMa.t[1];
		} else {
			rp = ml.aSPMb.t[0];
			rs = ml.aSPMa.t[0];
		}
		line = 0;
		if (Up) {
			if (lay == 0)
				q = i_ * ml.ki.t[lay] * (2. * ml.Zi.t[lay] - zj + zi)
						* gaussLeg.t[k][0];
			else
				q = i_ * ml.ki.t[lay] * fabs(zj - zi) * gaussLeg.t[k][0];
		} else {
			if (lay == ml.nbrMediums - 1)
				q = i_ * ml.ki.t[lay] * (zj + zi - 2. * ml.Zi.t[lay - 1])
						* gaussLeg.t[k][0];
			else
				q = i_ * ml.ki.t[lay] * fabs(zj - zi) * gaussLeg.t[k][0];
		}
		rauij = sqrt(dx * dx + dy * dy);
		eiphi = (dx + i_ * dy) / rauij;
		if (rauij == 0.)
			eiphi = 1.;
		bes.assign(N,
				ml.ki.t[lay] * rauij
						* sqrt(1. - gaussLeg.t[k][0] * gaussLeg.t[k][0]));
		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				col = 0;
				for (int l = 1; l <= N; l++) {
					for (int kk = -l; kk <= l; kk++) {
						if (Up) {
							mult = gaussLeg.t[k][1] * 2. * powm1(kk - l)
									* exp(q) * powi(l - n + abs(m - kk))
									* pow(eiphi, double(m - kk))
									* bes.Jn(abs(m - kk))
									/ sqrt(
											double(n * n + n)
													* double(l * l + l));

							RR.t[line][col] -= mult
									* (double(m * kk) * leg[k].Pimn(n, m)
											* leg[k].Pimn(l, kk) * rp
											- leg[k].Taumn(n, m)
													* leg[k].Taumn(l, kk) * rs);
							RR.t[line][col + NM] -= mult
									* (-double(m) * leg[k].Pimn(n, m)
											* leg[k].Taumn(l, kk) * rp
											+ double(kk) * leg[k].Taumn(n, m)
													* leg[k].Pimn(l, kk) * rs);
							RR.t[line + NM][col] -= mult
									* (double(kk) * leg[k].Taumn(n, m)
											* leg[k].Pimn(l, kk) * rp
											- double(m) * leg[k].Pimn(n, m)
													* leg[k].Taumn(l, kk) * rs);
							RR.t[line + NM][col + NM] -= mult
									* (-leg[k].Taumn(n, m) * leg[k].Taumn(l, kk)
											* rp
											+ double(m * kk) * leg[k].Pimn(n, m)
													* leg[k].Pimn(l, kk) * rs);
						} else {
							mult = gaussLeg.t[k][1] * 2. * powm1(m - n) * exp(q)
									* powi(l - n + abs(m - kk))
									* pow(eiphi, double(m - kk))
									* bes.Jn(abs(m - kk))
									/ sqrt(
											double(n * n + n)
													* double(l * l + l));
							RR.t[line][col] -= mult
									* (double(m * kk) * leg[k].Pimn(n, m)
											* leg[k].Pimn(l, kk) * rp
											- leg[k].Taumn(n, m)
													* leg[k].Taumn(l, kk) * rs);
							RR.t[line][col + NM] -= mult
									* (double(m) * leg[k].Pimn(n, m)
											* leg[k].Taumn(l, kk) * rp
											- double(kk) * leg[k].Taumn(n, m)
													* leg[k].Pimn(l, kk) * rs);
							RR.t[line + NM][col] -= mult
									* (-double(kk) * leg[k].Taumn(n, m)
											* leg[k].Pimn(l, kk) * double(m)
											* rp
											+ leg[k].Pimn(n, m)
													* leg[k].Taumn(l, kk) * rs);
							RR.t[line + NM][col + NM] -= mult
									* (-leg[k].Taumn(n, m) * leg[k].Taumn(l, kk)
											* rp
											+ double(m * kk) * leg[k].Pimn(n, m)
													* leg[k].Pimn(l, kk) * rs);
						}

						col++;
					}
				}
				line++;
			}
		}
	}

	for (int k = 0; k < Ng; k++) {
//    ml.setSourceInt(1.,0.,1.,0.,acos(xInt[k]),0.,zi);
		if (Up) {
			rp = ml.aSPMb.t[1];
			rs = ml.aSPMa.t[1];
		} else {
			rp = ml.aSPMb.t[0];
			rs = ml.aSPMa.t[0];
		}
		line = 0;
		if (Up) {
			if (lay == 0)
				q = i_ * ml.ki.t[lay] * fabs(2. * ml.Zi.t[lay] - zj + zi)
						* xInt[k];
			else
				q = i_ * ml.ki.t[lay] * fabs(zj - zi) * xInt[k];
		} else {
			if (lay == ml.nbrMediums - 1)
				q = i_ * ml.ki.t[lay] * fabs(zj + zi - 2. * ml.Zi.t[lay - 1])
						* xInt[k];
			else
				q = i_ * ml.ki.t[lay] * fabs(zj - zi) * xInt[k];
		}

		rauij = sqrt(dx * dx + dy * dy);
		eiphi = (dx + i_ * dy) / rauij;
		if (rauij == 0.)
			eiphi = 1.;
		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				col = 0;
				for (int l = 1; l <= N; l++) {
					for (int kk = -l; kk <= l; kk++) {

						if (Up) {
							mult = pasX * 2. * i_ * powm1(kk - l) * exp(q)
									* powi(l - n + abs(m - kk))
									* pow(eiphi, double(m - kk))
									* bes.Jn(abs(m - kk))
									/ sqrt(
											double(n * n + n)
													* double(l * l + l));
							if ((k == 0) || (k == Ng - 1))
								mult *= 0.5;
							RR.t[line][col] +=
									mult
											* (double(m * kk)
													* legC[k].Pimn(n, m)
													* legC[k].Pimn(l, kk) * rp
													- legC[k].Taumn(n, m)
															* legC[k].Taumn(l,
																	kk) * rs);
							RR.t[line][col + NM] += mult
									* (-double(m) * legC[k].Pimn(n, m)
											* legC[k].Taumn(l, kk) * rp
											+ double(kk) * legC[k].Taumn(n, m)
													* legC[k].Pimn(l, kk) * rs);
							RR.t[line + NM][col] +=
									mult
											* (double(kk) * legC[k].Taumn(n, m)
													* legC[k].Pimn(l, kk) * rp
													- double(m)
															* legC[k].Pimn(n, m)
															* legC[k].Taumn(l,
																	kk) * rs);
							RR.t[line + NM][col + NM] += mult
									* (-legC[k].Taumn(n, m)
											* legC[k].Taumn(l, kk) * rp
											+ double(m * kk)
													* legC[k].Pimn(n, m)
													* legC[k].Pimn(l, kk) * rs);
						} else {
							mult = pasX * 2. * i_ * powm1(m - n) * exp(q)
									* powi(l - n + abs(m - kk))
									* pow(eiphi, double(m - kk))
									* bes.Jn(abs(m - kk))
									/ sqrt(
											double(n * n + n)
													* double(l * l + l));
							if ((k == 0) || (k == Ng - 1))
								mult *= 0.5;
							RR.t[line][col] +=
									mult
											* (double(m * kk)
													* legC[k].Pimn(n, m)
													* legC[k].Pimn(l, kk) * rp
													- legC[k].Taumn(n, m)
															* legC[k].Taumn(l,
																	kk) * rs);
							RR.t[line][col + NM] += mult
									* (double(m) * legC[k].Pimn(n, m)
											* legC[k].Taumn(l, kk) * rp
											- double(kk) * legC[k].Taumn(n, m)
													* legC[k].Pimn(l, kk) * rs);
							RR.t[line + NM][col] +=
									mult
											* (-double(kk) * legC[k].Taumn(n, m)
													* legC[k].Pimn(l, kk)
													* double(m) * rp
													+ legC[k].Pimn(n, m)
															* legC[k].Taumn(l,
																	kk) * rs);
							RR.t[line + NM][col + NM] += mult
									* (-legC[k].Taumn(n, m)
											* legC[k].Taumn(l, kk) * rp
											+ double(m * kk)
													* legC[k].Pimn(n, m)
													* legC[k].Pimn(l, kk) * rs);
						}
						col++;
					}
				}
				line++;
			}
		}
	}

	return RR;
}
