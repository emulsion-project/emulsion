#include "AxisSymParticle.h"
#include <iostream>
#include <fstream>
using namespace std;

AxisSymParticle::AxisSymParticle() {
	T = R = NULL;
}

void AxisSymParticle::setParticleProperties(const int &typ, const int &mat,
		const int &matExt, const double &lambd, const RVector &rr) {
	type = typ;
	material = mat;
	materialExt = matExt;
	lambda = lambd;
	geom = rr;
	RefractiveIndex ri(material, 50e-9);
	RefractiveIndex re(matExt, 1e-9);
	n2 = ri.calcIndex(lambda);
	k2 = 2. * piDbl * n2 / lambda;
	n1 = re.calcIndex(lambda);
	k1 = 2. * piDbl * n1 / lambda;
}

void AxisSymParticle::setLambda(const double &lambd) {
	lambda = lambd;
	RefractiveIndex ri(material, 50e-9);
	RefractiveIndex re(materialExt, 1e-9);
	n2 = ri.calcIndex(lambda);
	k2 = 2. * piDbl * n2 / lambda;
	n1 = re.calcIndex(lambda);
	k1 = 2. * piDbl * n1 / lambda;
}

double AxisSymParticle::calcR(const double &cost) {
	double thet;
	if (type == 0)
		return geom.t[0];
	else if (type == 1)
		return geom.t[1]
				/ sqrt(
						1.
								+ (geom.t[1] * geom.t[1]
										/ (geom.t[0] * geom.t[0]) - 1.) * cost
										* cost);
	else if (type == 2) {
		thet = acos(cost);
		if (thet < atan(geom.t[1] / geom.t[0]))
			return geom.t[0] / cost;
		else if (thet < piDbl - atan(geom.t[1] / geom.t[0]))
			return geom.t[1] / sin(thet);
		else
			return -geom.t[0] / cost;
	}
	return 0.;
}

double AxisSymParticle::calcDR(const double &cost) {
	double thet;
	if (type == 0)
		return 0.;
	else if (type == 1) {
		if (geom.t[0] == geom.t[1])
			return 0.;
		else
			return geom.t[1]
					* (geom.t[1] * geom.t[1] / (geom.t[0] * geom.t[0]) - 1.)
					* cost * sqrt(1. - cost * cost)
					/ pow(
							1.
									+ (geom.t[1] * geom.t[1]
											/ (geom.t[0] * geom.t[0]) - 1.)
											* cost * cost, 1.5);
	} else if (type == 2) {
		thet = acos(cost);
		if (thet < atan(geom.t[1] / geom.t[0]))
			return geom.t[0] * sin(thet) / (cost * cost);
		else if (thet < piDbl - atan(geom.t[1] / geom.t[0]))
			return -geom.t[0] * cost / (sin(thet) * sin(thet));
		else
			return -geom.t[0] * sin(thet) / (cost * cost);
	}
	return 0.;
}

void AxisSymParticle::setTMatrixParameters(const int &n, const int &ng,
		const int &meth) {
	N = n;
	Ng = ng;
	method = meth;
	if (T != NULL)
		delete[] T;
	if (R != NULL)
		delete[] R;
	T = new Matrix[2 * N + 1];
	R = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		T[i] = Matrix(2 * nbr, 2 * nbr, 0.);
		R[i] = Matrix(2 * nbr, 2 * nbr, 0.);
	}
	if (type == 0)
		calcTSphere();
	else {
		if (method == 0) {
			if (type == 2)
				calcTLScylinder();
			else
				calcTEBCM();
		}
		if (method == 1)
			calcTBC();
		if (method == 2)
			calcTLS();
	}
	return;
}

void AxisSymParticle::calcTSphere() {
	int l, i = 0;
	Complex m2 = n2 / n1;
	SphericalBessel bes1(N, k1 * geom.t[0]);
	SphericalBessel bes2(N, k2 * geom.t[0]);
	for (int m = -N; m <= N; m++) {
		l = 0;
		nbr = N - abs(m) + 1;
		if (abs(m) == 0)
			nbr--;
		for (int n = abs(m); n <= N; n++) {
			if (n == 0)
				n++;
			T[i].t[l + nbr][l + nbr] = -(m2 * m2 * bes2.jn(n) * bes1.djn(n)
					- bes2.djn(n) * bes1.jn(n))
					/ (m2 * m2 * bes2.jn(n) * bes1.dh1n(n)
							- bes2.djn(n) * bes1.h1n(n)); // -an
			T[i].t[l][l] =
					-(bes1.jn(n) * bes2.djn(n) - bes1.djn(n) * bes2.jn(n))
							/ (bes2.djn(n) * bes1.h1n(n)
									- bes2.jn(n) * bes1.dh1n(n)); // -bn

			//R[i].t[l][l] = -(bes1.dh1n(n)*bes1.jn(n)-bes1.h1n(n)*bes1.djn(n))/(bes2.jn(n)*bes1.dh1n(n)-bes2.djn(n)*bes1.h1n(n)); // -cn
			//R[i].t[l+nbr][l+nbr] = -(m2*bes1.jn(n)*bes1.dh1n(n)-m2*bes1.djn(n)*bes1.h1n(n))/(m2*m2*bes2.jn(n)*bes1.dh1n(n)-bes2.djn(n)*bes1.h1n(n)/m2); // -dn

			R[i].t[l][l] = -i_ * m2
					/ (k2 * geom.t[0] * bes2.jn(n) * bes1.djn(n)
							- m2 * bes2.djn(n) * k1 * geom.t[0] * bes1.jn(n)); // -cn
			R[i].t[l + nbr][l + nbr] = -i_ * m2
					/ (m2 * k2 * geom.t[0] * bes2.jn(n) * bes1.djn(n)
							- bes2.djn(n) * k1 * geom.t[0] * bes1.jn(n)); // -dn
			l++;
		}
		i++;
	}
	return;
}

void AxisSymParticle::calcTEBCM() {
	int nbr;
	Matrix *Q11, *Q31;
	Q11 = new Matrix[2 * N + 1];
	Q31 = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		Q11[i] = Matrix(2 * nbr, 2 * nbr, 0.);
		Q31[i] = Matrix(2 * nbr, 2 * nbr, 0.);
	}
	RMatrix gaussLeg = gaussLegQuad(Ng, -1., 1.);
	int l, c;
	Vector3d ndS;
	VSWF vswf1, vswf2;
	double r, dr, cost;
	Complex w;
	for (int k = 0; k < Ng; k++) {
		cost = gaussLeg.t[k][0];
		w = gaussLeg.t[k][1] * i_ * k1 * k1 * 2.;
		r = calcR(cost);
		dr = calcDR(cost);
		ndS.t[0] = r * r;
		ndS.t[1] = -r * dr;
		ndS.t[2] = 0.;
		vswf1.setParameters(N, k1 * r, cost, 0.);
		vswf2.setParameters(N, k2 * r, cost, 0.);
		for (int m = -N; m <= N; m++) {
			l = 0;
			nbr = N - abs(m) + 1;
			if (abs(m) == 0)
				nbr--;
			for (int nl = abs(m); nl <= N; nl++) {
				if (nl == 0)
					nl++;
				c = 0;
				for (int nc = abs(m); nc <= N; nc++) {
					if (nc == 0)
						nc++;
					vswf1.calcVSWFs(-m, nl, 1, false);
					vswf2.calcVSWFs(m, nc, 1, false);
					Q11[m + N].t[l][c] += w
							* (dotProduct(crossProduct(vswf2.Mmn, vswf1.Nmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Nmn,
															vswf1.Mmn), ndS));
					Q11[m + N].t[l][c + nbr] += w
							* (dotProduct(crossProduct(vswf2.Nmn, vswf1.Nmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Mmn,
															vswf1.Mmn), ndS));
					Q11[m + N].t[l + nbr][c] += w
							* (dotProduct(crossProduct(vswf2.Mmn, vswf1.Mmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Nmn,
															vswf1.Nmn), ndS));
					Q11[m + N].t[l + nbr][c + nbr] += w
							* (dotProduct(crossProduct(vswf2.Nmn, vswf1.Mmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Mmn,
															vswf1.Nmn), ndS));

					vswf1.calcVSWFs(-m, nl, 3, false);
					Q31[m + N].t[l][c] -= w
							* (dotProduct(crossProduct(vswf2.Mmn, vswf1.Nmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Nmn,
															vswf1.Mmn), ndS));
					Q31[m + N].t[l][c + nbr] -= w
							* (dotProduct(crossProduct(vswf2.Nmn, vswf1.Nmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Mmn,
															vswf1.Mmn), ndS));
					Q31[m + N].t[l + nbr][c] -= w
							* (dotProduct(crossProduct(vswf2.Mmn, vswf1.Mmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Nmn,
															vswf1.Nmn), ndS));
					Q31[m + N].t[l + nbr][c + nbr] -= w
							* (dotProduct(crossProduct(vswf2.Nmn, vswf1.Mmn),
									ndS)
									+ n2 / n1
											* dotProduct(
													crossProduct(vswf2.Mmn,
															vswf1.Nmn), ndS));
					c++;
				}
				l++;
			}
		}
	}
	for (int i = 0; i < 2 * N + 1; i++) {
		invGaussj(Q31[i]);
		T[i] = Q11[i] * Q31[i];
		R[i] = Q31[i];
	}
	delete[] Q11;
	delete[] Q31;
	return;
}

void AxisSymParticle::calcTBC() {
	int nbr;
	Matrix *A, *B, *TBC;
	A = new Matrix[2 * N + 1];
	B = new Matrix[2 * N + 1];
	TBC = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		A[i] = Matrix(4 * nbr, 4 * nbr, 0.);
		B[i] = Matrix(4 * nbr, 2 * nbr, 0.);
		TBC[i] = Matrix(4 * nbr, 2 * nbr, 0.);
	}
	RMatrix gaussLeg = gaussLegQuad(Ng, -1., 1.);
	int l, c;
	Vector3d ndS;
	VSWF vswf1, vswf2, vswf2c;
	double r, dr, cost;
	Complex w;
	for (int k = 0; k < Ng; k++) {

		cost = gaussLeg.t[k][0];
		w = gaussLeg.t[k][1];
		r = calcR(cost);
		dr = calcDR(cost);
		ndS.t[0] = r * r;
		ndS.t[1] = -r * dr;
		ndS.t[2] = 0.;
		vswf1.setParameters(N, k1 * r, cost, 0.);
		vswf2.setParameters(N, k2 * r, cost, 0.);
		vswf2c.setParameters(N, k2 * r, cost, 0.);
		for (int m = -N; m <= N; m++) {
			l = 0;
			nbr = N - abs(m) + 1;
			if (abs(m) == 0)
				nbr--;
			for (int nl = abs(m); nl <= N; nl++) {
				if (nl == 0)
					nl++;
				c = 0;
				for (int nc = abs(m); nc <= N; nc++) {
					if (nc == 0)
						nc++;
					vswf1.calcVSWFs(m, nc, 3, false);
					vswf2.calcVSWFs(-m, nl, 1, false);
					A[m + N].t[l][c] += w * k2
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l][c + nbr] += w * k2
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l + nbr][c] += w * k2
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + nbr][c + nbr] += w * k2
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Mmn),
									ndS);

					A[m + N].t[l + 2 * nbr][c] += w * k1
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + 2 * nbr][c + nbr] += w * k1
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + 3 * nbr][c] += w * k1
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l + 3 * nbr][c + nbr] += w * k1
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Nmn),
									ndS);

					vswf2c.calcVSWFs(m, nc, 1, false);

					A[m + N].t[l][c + 2 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Mmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l][c + 3 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Nmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l + nbr][c + 2 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Mmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + nbr][c + 3 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Nmn, vswf2.Mmn),
									ndS);

					A[m + N].t[l + 2 * nbr][c + 2 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Nmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + 2 * nbr][c + 3 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Mmn, vswf2.Mmn),
									ndS);
					A[m + N].t[l + 3 * nbr][c + 2 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Nmn, vswf2.Nmn),
									ndS);
					A[m + N].t[l + 3 * nbr][c + 3 * nbr] -= w * k2
							* dotProduct(crossProduct(vswf2c.Mmn, vswf2.Nmn),
									ndS);

					vswf1.calcVSWFs(m, nc, 1, false);

					B[m + N].t[l][c] -= w * k2
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Nmn),
									ndS);
					B[m + N].t[l][c + nbr] -= w * k2
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Nmn),
									ndS);
					B[m + N].t[l + nbr][c] -= w * k2
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Mmn),
									ndS);
					B[m + N].t[l + nbr][c + nbr] -= w * k2
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Mmn),
									ndS);

					B[m + N].t[l + 2 * nbr][c] -= w * k1
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Mmn),
									ndS);
					B[m + N].t[l + 2 * nbr][c + nbr] -= w * k1
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Mmn),
									ndS);
					B[m + N].t[l + 3 * nbr][c] -= w * k1
							* dotProduct(crossProduct(vswf1.Nmn, vswf2.Nmn),
									ndS);
					B[m + N].t[l + 3 * nbr][c + nbr] -= w * k1
							* dotProduct(crossProduct(vswf1.Mmn, vswf2.Nmn),
									ndS);
					c++;
				}
				l++;
			}
		}
	}

	for (int i = 0; i < 2 * N + 1; i++) {
		invGaussj(A[i]);
		TBC[i] = A[i] * B[i];
		for (int n = 0; n < TBC[i].m; n++) {
			for (int m = 0; m < TBC[i].m; m++) {
				T[i].t[n][m] = TBC[i].t[n][m];
				R[i].t[n][m] = TBC[i].t[n + TBC[i].m][m];
			}
		}
	}
	delete[] A;
	delete[] B;
	delete[] TBC;
	return;
}

/*void AxisSymParticle::calcTLS() {
 int nbr;
 Matrix *A,*B,*TBC;
 A = new Matrix[2*N+1];
 B = new Matrix[2*N+1];
 TBC = new Matrix[2*N+1];
 for (int i=0 ; i<2*N+1 ; i++) {
 nbr = N-abs(i-N)+1;
 if (abs(i-N)==0) nbr--;
 A[i] = Matrix(4*nbr,4*nbr,0.);
 B[i] = Matrix(4*nbr,2*nbr,0.);
 TBC[i] = Matrix(4*nbr,2*nbr,0.);
 }
 RMatrix gaussLeg = gaussLegQuad(Ng,-1.,1.);
 int l,c;
 Vector3d ndS;
 VSWF vswf1,vswf1c,vswf2,vswf2c;
 double r,dr,cost,dS;
 Complex w;
 for (int k=0 ; k<Ng ; k++) {

 cost = gaussLeg.t[k][0];
 w = gaussLeg.t[k][1];
 r = calcR(cost);
 dr = calcDR(cost);
 ndS.t[0] = r*r;
 ndS.t[1] = -r*dr;
 ndS.t[2] = 0.;

 vswf1.setParameters(N,k1*r,cost,0.);
 vswf1c.setParameters(N,k1*r,cost,0.);
 vswf2.setParameters(N,k2*r,cost,0.);
 vswf2c.setParameters(N,k2*r,cost,0.);
 for (int m=-N ; m<=N ; m++) {
 l = 0;
 nbr = N-abs(m)+1;
 if (abs(m)==0) nbr--;
 for (int nl=abs(m) ; nl<=N ; nl++) {
 if (nl==0) nl++;
 c = 0;
 for (int nc=abs(m) ; nc<=N ; nc++) {
 if (nc==0) nc++;
 vswf1.calcVSWFs(m,nc,3,false);
 vswf2.calcVSWFs(-m,nl,1,false);
 //A[m+N].t[l][c] += w*k2*dotProduct(crossProduct(vswf1.Mmn,vswf2.Nmn),ndS);
 //A[m+N].t[l][c+nbr] += w*k2*dotProduct(crossProduct(vswf1.Nmn,vswf2.Nmn),ndS);
 //A[m+N].t[l+nbr][c] += w*k2*dotProduct(crossProduct(vswf1.Mmn,vswf2.Mmn),ndS);
 //A[m+N].t[l+nbr][c+nbr] += w*k2*dotProduct(crossProduct(vswf1.Nmn,vswf2.Mmn),ndS);

 A[m+N].t[l][c] += w*k1*dotProduct(crossProduct(vswf1.Nmn,vswf2.Mmn),ndS);
 A[m+N].t[l][c+nbr] += w*k1*dotProduct(crossProduct(vswf1.Mmn,vswf2.Mmn),ndS);
 A[m+N].t[l+nbr][c] += w*k1*dotProduct(crossProduct(vswf1.Nmn,vswf2.Nmn),ndS);
 A[m+N].t[l+nbr][c+nbr] += w*k1*dotProduct(crossProduct(vswf1.Mmn,vswf2.Nmn),ndS);

 vswf2c.calcVSWFs(m,nc,1,false);

 //A[m+N].t[l][c+2*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Mmn,vswf2.Nmn),ndS);
 //A[m+N].t[l][c+3*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Nmn,vswf2.Nmn),ndS);
 //A[m+N].t[l+nbr][c+2*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Mmn,vswf2.Mmn),ndS);
 //A[m+N].t[l+nbr][c+3*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Nmn,vswf2.Mmn),ndS);

 A[m+N].t[l][c+2*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Nmn,vswf2.Mmn),ndS);
 A[m+N].t[l][c+3*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Mmn,vswf2.Mmn),ndS);
 A[m+N].t[l+nbr][c+2*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Nmn,vswf2.Nmn),ndS);
 A[m+N].t[l+nbr][c+3*nbr] -= w*k2*dotProduct(crossProduct(vswf2c.Mmn,vswf2.Nmn),ndS);

 vswf1.calcVSWFs(m,nc,1,false);

 //B[m+N].t[l][c] -= w*k2*dotProduct(crossProduct(vswf1.Mmn,vswf2.Nmn),ndS);
 //B[m+N].t[l][c+nbr] -= w*k2*dotProduct(crossProduct(vswf1.Nmn,vswf2.Nmn),ndS);
 //B[m+N].t[l+nbr][c] -= w*k2*dotProduct(crossProduct(vswf1.Mmn,vswf2.Mmn),ndS);
 //B[m+N].t[l+nbr][c+nbr] -= w*k2*dotProduct(crossProduct(vswf1.Nmn,vswf2.Mmn),ndS);

 B[m+N].t[l][c] -= w*k1*dotProduct(crossProduct(vswf1.Nmn,vswf2.Mmn),ndS);
 B[m+N].t[l][c+nbr] -= w*k1*dotProduct(crossProduct(vswf1.Mmn,vswf2.Mmn),ndS);
 B[m+N].t[l+nbr][c] -= w*k1*dotProduct(crossProduct(vswf1.Nmn,vswf2.Nmn),ndS);
 B[m+N].t[l+nbr][c+nbr] -= w*k1*dotProduct(crossProduct(vswf1.Mmn,vswf2.Nmn),ndS);


 vswf1c.calcVSWFs(m,nc,3,false);
 vswf1.calcVSWFs(m,nl,3,false);
 vswf2.calcVSWFs(m,nl,1,false);
 vswf2c.calcVSWFs(m,nc,1,false);
 A[m+N].t[l+2*nbr][c] -= w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Mmn)) + dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Nmn)));
 A[m+N].t[l+2*nbr][c+nbr] -= w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Mmn)) + dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Nmn)));
 A[m+N].t[l+3*nbr][c] -= w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Nmn)) + dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Mmn)));
 A[m+N].t[l+3*nbr][c+nbr] -= w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Nmn)) + dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Mmn)));

 //A[m+N].t[l+2*nbr][c] += w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Nmn)));
 //A[m+N].t[l+2*nbr][c+nbr] += w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Nmn)));
 //A[m+N].t[l+3*nbr][c] += w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Mmn)));
 //A[m+N].t[l+3*nbr][c+nbr] += w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Mmn)));

 A[m+N].t[l+2*nbr][c+2*nbr] += w*(dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf1.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf1.Nmn)));
 A[m+N].t[l+2*nbr][c+3*nbr] += w*(dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf1.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf1.Nmn)));
 A[m+N].t[l+3*nbr][c+2*nbr] += w*(dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf1.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf1.Mmn)));
 A[m+N].t[l+3*nbr][c+3*nbr] += w*(dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf1.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf1.Mmn)));

 //A[m+N].t[l+2*nbr][c+2*nbr] -= w*(dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf2.Mmn)) + k2*k2/(k1*k1)*dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf2.Nmn)));
 //A[m+N].t[l+2*nbr][c+3*nbr] -= w*(dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf2.Mmn)) + k2*k2/(k1*k1)*dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf2.Nmn)));
 //A[m+N].t[l+3*nbr][c+2*nbr] -= w*(dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf2.Nmn)) + k2*k2/(k1*k1)*dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf2.Mmn)));
 //A[m+N].t[l+3*nbr][c+3*nbr] -= w*(dotProduct(crossProduct(ndS,vswf2c.Nmn),crossProduct(ndS,vswf2.Nmn)) + k2*k2/(k1*k1)*dotProduct(crossProduct(ndS,vswf2c.Mmn),crossProduct(ndS,vswf2.Mmn)));

 vswf1c.calcVSWFs(m,nc,1,false);

 B[m+N].t[l+2*nbr][c] += w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Mmn)) + dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Nmn)));
 B[m+N].t[l+2*nbr][c+nbr] += w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Mmn)) + dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Nmn)));
 B[m+N].t[l+3*nbr][c] += w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Nmn)) + dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Mmn)));
 B[m+N].t[l+3*nbr][c+nbr] += w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf1.Nmn)) + dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf1.Mmn)));


 //B[m+N].t[l+2*nbr][c] -= w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Nmn)));
 //B[m+N].t[l+2*nbr][c+nbr] -= w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Mmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Nmn)));
 //B[m+N].t[l+3*nbr][c] -= w*(dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Mmn)));
 //B[m+N].t[l+3*nbr][c+nbr] -= w*(dotProduct(crossProduct(ndS,vswf1c.Nmn),crossProduct(ndS,vswf2.Nmn)) + k2/k1*dotProduct(crossProduct(ndS,vswf1c.Mmn),crossProduct(ndS,vswf2.Mmn)));
 c++;
 }
 l++;
 }
 }
 }
 //ofstream of("outMat.txt");
 for (int i=0 ; i<2*N+1 ; i++) {
 invGaussj(A[i]);
 //of << A[i] << endl;
 //TBC[i] = invCholesky(A[i])*B[i];
 TBC[i] = A[i]*B[i];
 for (int n=0 ; n<TBC[i].m ; n++) {
 for (int m=0 ; m<TBC[i].m ; m++) {
 T[i].t[n][m] = TBC[i].t[n][m];
 R[i].t[n][m] = TBC[i].t[n+TBC[i].m][m];
 }
 }
 }
 delete [] A;
 delete [] B;
 delete [] TBC;
 return;
 }*/

void AxisSymParticle::calcTLS() {
	int nbr;
	Matrix *A, *B, *TBC;
	A = new Matrix[2 * N + 1];
	B = new Matrix[2 * N + 1];
	TBC = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		A[i] = Matrix(4 * nbr, 4 * nbr, 0.);
		B[i] = Matrix(4 * nbr, 2 * nbr, 0.);
		TBC[i] = Matrix(4 * nbr, 2 * nbr, 0.);
	}
	RMatrix gaussLeg = gaussLegQuad(Ng, -1., 1.);
	int l, c;
	Vector3d ndS;
	VSWF vswf1, vswf1c, vswf2, vswf2c;
	double r, dr, cost, dS;
	Complex w;
	for (int k = 0; k < Ng; k++) {

		cost = gaussLeg.t[k][0];
		w = gaussLeg.t[k][1];
		r = calcR(cost);
		dr = calcDR(cost);
		ndS.t[0] = r * r;
		ndS.t[1] = -r * dr;
		ndS.t[2] = 0.;

		vswf1.setParameters(N, k1 * r, cost, 0.);
		vswf1c.setParameters(N, k1 * r, cost, 0.);
		vswf2.setParameters(N, k2 * r, cost, 0.);
		vswf2c.setParameters(N, k2 * r, cost, 0.);
		for (int m = -N; m <= N; m++) {
			l = 0;
			nbr = N - abs(m) + 1;
			if (abs(m) == 0)
				nbr--;
			for (int nl = abs(m); nl <= N; nl++) {
				if (nl == 0)
					nl++;
				c = 0;
				for (int nc = abs(m); nc <= N; nc++) {
					if (nc == 0)
						nc++;
					vswf1c.calcVSWFs(m, nc, 3, false);
					vswf1.calcVSWFs(m, nl, 3, false);
					vswf2.calcVSWFs(m, nl, 1, false);
					vswf2c.calcVSWFs(m, nc, 1, false);
					A[m + N].t[l][c] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf1.Mmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Nmn),
											crossProduct(ndS, vswf1.Nmn)));
					A[m + N].t[l][c + nbr] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf1.Mmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Mmn),
											crossProduct(ndS, vswf1.Nmn)));
					A[m + N].t[l + nbr][c] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf1.Nmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Nmn),
											crossProduct(ndS, vswf1.Mmn)));
					A[m + N].t[l + nbr][c + nbr] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf1.Nmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Mmn),
											crossProduct(ndS, vswf1.Mmn)));

					A[m + N].t[l + 2 * nbr][c] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Nmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					A[m + N].t[l + 2 * nbr][c + nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Mmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					A[m + N].t[l + 3 * nbr][c] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Nmn),
													crossProduct(ndS,
															vswf2.Mmn)));
					A[m + N].t[l + 3 * nbr][c + nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Mmn),
													crossProduct(ndS,
															vswf2.Mmn)));

					A[m + N].t[l][c + 2 * nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf2c.Mmn),
									crossProduct(ndS, vswf1.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Nmn),
													crossProduct(ndS,
															vswf1.Nmn)));
					A[m + N].t[l][c + 3 * nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf2c.Nmn),
									crossProduct(ndS, vswf1.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Mmn),
													crossProduct(ndS,
															vswf1.Nmn)));
					A[m + N].t[l + nbr][c + 2 * nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf2c.Mmn),
									crossProduct(ndS, vswf1.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Nmn),
													crossProduct(ndS,
															vswf1.Mmn)));
					A[m + N].t[l + nbr][c + 3 * nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf2c.Nmn),
									crossProduct(ndS, vswf1.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Mmn),
													crossProduct(ndS,
															vswf1.Mmn)));

					A[m + N].t[l + 2 * nbr][c + 2 * nbr] += w
							* (dotProduct(crossProduct(ndS, vswf2c.Mmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 * k2 / (k1 * k1)
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Nmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					A[m + N].t[l + 2 * nbr][c + 3 * nbr] += w
							* (dotProduct(crossProduct(ndS, vswf2c.Nmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 * k2 / (k1 * k1)
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Mmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					A[m + N].t[l + 3 * nbr][c + 2 * nbr] += w
							* (dotProduct(crossProduct(ndS, vswf2c.Mmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 * k2 / (k1 * k1)
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Nmn),
													crossProduct(ndS,
															vswf2.Mmn)));
					A[m + N].t[l + 3 * nbr][c + 3 * nbr] += w
							* (dotProduct(crossProduct(ndS, vswf2c.Nmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 * k2 / (k1 * k1)
											* dotProduct(
													crossProduct(ndS,
															vswf2c.Mmn),
													crossProduct(ndS,
															vswf2.Mmn)));

					vswf1c.calcVSWFs(m, nc, 1, false);

					B[m + N].t[l][c] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf1.Mmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Nmn),
											crossProduct(ndS, vswf1.Nmn)));
					B[m + N].t[l][c + nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf1.Mmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Mmn),
											crossProduct(ndS, vswf1.Nmn)));
					B[m + N].t[l + nbr][c] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf1.Nmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Nmn),
											crossProduct(ndS, vswf1.Mmn)));
					B[m + N].t[l + nbr][c + nbr] -= w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf1.Nmn))
									+ dotProduct(crossProduct(ndS, vswf1c.Mmn),
											crossProduct(ndS, vswf1.Mmn)));

					B[m + N].t[l + 2 * nbr][c] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Nmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					B[m + N].t[l + 2 * nbr][c + nbr] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf2.Mmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Mmn),
													crossProduct(ndS,
															vswf2.Nmn)));
					B[m + N].t[l + 3 * nbr][c] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Mmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Nmn),
													crossProduct(ndS,
															vswf2.Mmn)));
					B[m + N].t[l + 3 * nbr][c + nbr] += w
							* (dotProduct(crossProduct(ndS, vswf1c.Nmn),
									crossProduct(ndS, vswf2.Nmn))
									+ k2 / k1
											* dotProduct(
													crossProduct(ndS,
															vswf1c.Mmn),
													crossProduct(ndS,
															vswf2.Mmn)));
					c++;
				}
				l++;
			}
		}
	}
	ofstream of("outMat.txt");
	for (int i = 0; i < 2 * N + 1; i++) {
		invGaussj(A[i]);
		//of << A[i] << endl;
		//TBC[i] = invCholesky(A[i])*B[i];
		TBC[i] = A[i] * B[i];
		for (int n = 0; n < TBC[i].m; n++) {
			for (int m = 0; m < TBC[i].m; m++) {
				T[i].t[n][m] = TBC[i].t[n][m];
				R[i].t[n][m] = TBC[i].t[n + TBC[i].m][m];
			}
		}
	}
	delete[] A;
	delete[] B;
	delete[] TBC;
	return;
}

void AxisSymParticle::calcTLScylinder() {
	int nbr;
	Matrix *Q11, *Q31;
	Q11 = new Matrix[2 * N + 1];
	Q31 = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		Q11[i] = Matrix(2 * nbr, 2 * nbr, 0.);
		Q31[i] = Matrix(2 * nbr, 2 * nbr, 0.);
	}
	RMatrix gaussLeg;
	int l, c;
	Vector3d ndS;
	VSWF vswf1, vswf2;
	double r, dr, cost, w;
	for (int interv = 0; interv < 3; interv++) {
		if (interv == 0)
			gaussLeg = gaussLegQuad(Ng, -1.,
					cos(piDbl - atan(geom.t[1] / geom.t[0])));
		if (interv == 1)
			gaussLeg = gaussLegQuad(Ng,
					cos(piDbl - atan(geom.t[1] / geom.t[0])),
					cos(atan(geom.t[1] / geom.t[0])));
		if (interv == 2)
			gaussLeg = gaussLegQuad(Ng, cos(atan(geom.t[1] / geom.t[0])), 1.);
		for (int k = 0; k < Ng; k++) {
			cost = gaussLeg.t[k][0];
			w = gaussLeg.t[k][1];
			r = calcR(cost);

			dr = calcDR(cost);
			ndS.t[0] = r * r;
			ndS.t[1] = -r * dr;
			ndS.t[2] = 0.;
			vswf1.setParameters(N, k1 * r, cost, 0.);
			vswf2.setParameters(N, k2 * r, cost, 0.);
			for (int m = -N; m <= N; m++) {
				l = 0;
				nbr = N - abs(m) + 1;
				if (abs(m) == 0)
					nbr--;
				for (int nl = abs(m); nl <= N; nl++) {
					if (nl == 0)
						nl++;
					c = 0;
					for (int nc = abs(m); nc <= N; nc++) {
						if (nc == 0)
							nc++;
						vswf1.calcVSWFs(-m, nl, 1, false);
						vswf2.calcVSWFs(m, nc, 1, false);
						Q11[m + N].t[l][c] +=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Mmn,
																vswf1.Nmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Nmn,
																		vswf1.Mmn),
																ndS));
						Q11[m + N].t[l][c + nbr] +=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Nmn,
																vswf1.Nmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Mmn,
																		vswf1.Mmn),
																ndS));
						Q11[m + N].t[l + nbr][c] +=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Mmn,
																vswf1.Mmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Nmn,
																		vswf1.Nmn),
																ndS));
						Q11[m + N].t[l + nbr][c + nbr] +=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Nmn,
																vswf1.Mmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Mmn,
																		vswf1.Nmn),
																ndS));

						vswf1.calcVSWFs(-m, nl, 3, false);
						Q31[m + N].t[l][c] -=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Mmn,
																vswf1.Nmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Nmn,
																		vswf1.Mmn),
																ndS));
						Q31[m + N].t[l][c + nbr] -=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Nmn,
																vswf1.Nmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Mmn,
																		vswf1.Mmn),
																ndS));
						Q31[m + N].t[l + nbr][c] -=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Mmn,
																vswf1.Mmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Nmn,
																		vswf1.Nmn),
																ndS));
						Q31[m + N].t[l + nbr][c + nbr] -=
								w
										* (k1
												* dotProduct(
														crossProduct(vswf2.Nmn,
																vswf1.Mmn), ndS)
												+ k2
														* dotProduct(
																crossProduct(
																		vswf2.Mmn,
																		vswf1.Nmn),
																ndS));
						c++;
					}
					l++;
				}
			}
		}
	}
	for (int i = 0; i < 2 * N + 1; i++) {
		invGaussj(Q31[i]);
		R[i] = Q31[i];
		T[i] = Q11[i] * Q31[i];
	}
	delete[] Q11;
	delete[] Q31;
	return;
}

void AxisSymParticle::calcSources() {
	sources = Vector(N, 0.);
	if (geom.t[0] >= geom.t[1]) {
		double d = 0.95 * geom.t[0];
		double step = 2. * d / double(N - 1);
		for (int i = 0; i < N; i++)
			sources.t[i] = (-d + double(i) * step);
	} else {
		double d = 2. * geom.t[1];
		double step = d / double(N + 1);
		for (int i = 0; i < N; i++)
			sources.t[i] = i_ * (-d * 0.5 + double(i + 1) * step);
	}
	return;
}

void AxisSymParticle::calcTDS() {
	int nbr;
	Matrix *Q11, *Q311, *Q312;
	Q11 = new Matrix[2 * N + 1];
	Q311 = new Matrix[2 * N + 1];
	Q312 = new Matrix[2 * N + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		nbr = N - abs(i - N) + 1;
		if (abs(i - N) == 0)
			nbr--;
		Q11[i] = Matrix(2 * nbr, 2 * N, 0.);
		Q311[i] = Matrix(2 * N, 2 * nbr, 0.);
		Q312[i] = Matrix(2 * N, 2 * N, 0.);
	}
	RMatrix gaussLeg = gaussLegQuad(Ng, 0., piDbl);
	int l, c;
	Vector3d ndS;
	VSWF vswf1;
	VSWFds vswf2ds, vswf1ds;
	Complex mult, rap = n2 / n1;
	double r, dr, theta, cost;
	Complex w;
	calcSources();
	for (int k = 0; k < Ng; k++) {
		theta = gaussLeg.t[k][0];
		cost = cos(theta);
		w = gaussLeg.t[k][1];
		mult = w * 2. * i_ * k1;
		r = calcR(cost);
		dr = calcDR(cost);
		ndS.t[0] = r * r * sin(theta);
		ndS.t[1] = -r * dr * sin(theta);
		ndS.t[2] = 0.;
		vswf2ds.setParameters(N, k2, r, sources, theta, 0.);
		vswf1ds.setParameters(N, k1, r, sources, theta, 0.);
		vswf1.setParameters(N, k1 * r, cost, 0.);
		for (int m = -N; m <= N; m++) {
			l = 0;
			nbr = N - abs(m) + 1;
			if (abs(m) == 0)
				nbr--;
			for (int nl = abs(m); nl <= N; nl++) {
				if (nl == 0)
					nl++;
				c = 0;
				for (int nc = 1; nc <= N; nc++) {
					vswf1.calcVSWFs(-m, nl, 1, false);
					vswf2ds.calcVSWFsds(m, nc, 1, false);
					Q11[m + N].t[l][c] += mult
							* (k1
									* dotProduct(
											crossProduct(vswf2ds.Mmn,
													vswf1.Nmn), ndS)
									+ k2
											* dotProduct(
													crossProduct(vswf2ds.Nmn,
															vswf1.Mmn), ndS));
					Q11[m + N].t[l][c + N] += mult
							* (k1
									* dotProduct(
											crossProduct(vswf2ds.Nmn,
													vswf1.Nmn), ndS)
									+ k2
											* dotProduct(
													crossProduct(vswf2ds.Mmn,
															vswf1.Mmn), ndS));
					Q11[m + N].t[l + nbr][c] += mult
							* (k1
									* dotProduct(
											crossProduct(vswf2ds.Mmn,
													vswf1.Mmn), ndS)
									+ k2
											* dotProduct(
													crossProduct(vswf2ds.Nmn,
															vswf1.Nmn), ndS));
					Q11[m + N].t[l + nbr][c + N] += mult
							* (k1
									* dotProduct(
											crossProduct(vswf2ds.Nmn,
													vswf1.Mmn), ndS)
									+ k2
											* dotProduct(
													crossProduct(vswf2ds.Mmn,
															vswf1.Nmn), ndS));
					c++;
				}
				l++;
			}
		}
		w *= k1 * k1;
		for (int m = -N; m <= N; m++) {
			l = 0;
			nbr = N - abs(m) + 1;
			if (abs(m) == 0)
				nbr--;
			for (int nl = 1; nl <= N; nl++) {
				c = 0;
				for (int nc = abs(m); nc <= N; nc++) {
					if (nc == 0)
						nc++;
					vswf1.calcVSWFs(m, nc, 1, false);
					vswf1ds.calcVSWFsds(-m, nl, 3, false);
					Q311[m + N].t[l][c] += w
							* (dotProduct(crossProduct(vswf1.Mmn, vswf1ds.Nmn),
									ndS)
									+ dotProduct(
											crossProduct(vswf1.Nmn,
													vswf1ds.Mmn), ndS));
					Q311[m + N].t[l][c + nbr] += w
							* (dotProduct(crossProduct(vswf1.Nmn, vswf1ds.Nmn),
									ndS)
									+ dotProduct(
											crossProduct(vswf1.Mmn,
													vswf1ds.Mmn), ndS));
					Q311[m + N].t[l + N][c] += w
							* (dotProduct(crossProduct(vswf1.Mmn, vswf1ds.Mmn),
									ndS)
									+ dotProduct(
											crossProduct(vswf1.Nmn,
													vswf1ds.Nmn), ndS));
					Q311[m + N].t[l + N][c + nbr] += w
							* (dotProduct(crossProduct(vswf1.Nmn, vswf1ds.Mmn),
									ndS)
									+ dotProduct(
											crossProduct(vswf1.Mmn,
													vswf1ds.Nmn), ndS));
					c++;
				}
				l++;
			}
		}
		for (int m = -N; m <= N; m++) {
			l = 0;
			for (int nl = 1; nl <= N; nl++) {
				c = 0;
				for (int nc = 1; nc <= N; nc++) {
					vswf2ds.calcVSWFsds(m, nc, 1, false);
					vswf1ds.calcVSWFsds(-m, nl, 3, false);
					Q312[m + N].t[l][c] += w
							* (dotProduct(
									crossProduct(vswf2ds.Mmn, vswf1ds.Nmn), ndS)
									+ rap
											* dotProduct(
													crossProduct(vswf2ds.Nmn,
															vswf1ds.Mmn), ndS));
					Q312[m + N].t[l][c + N] += w
							* (dotProduct(
									crossProduct(vswf2ds.Nmn, vswf1ds.Nmn), ndS)
									+ rap
											* dotProduct(
													crossProduct(vswf2ds.Mmn,
															vswf1ds.Mmn), ndS));
					Q312[m + N].t[l + N][c] += w
							* (dotProduct(
									crossProduct(vswf2ds.Mmn, vswf1ds.Mmn), ndS)
									+ rap
											* dotProduct(
													crossProduct(vswf2ds.Nmn,
															vswf1ds.Nmn), ndS));
					Q312[m + N].t[l + N][c + N] += w
							* (dotProduct(
									crossProduct(vswf2ds.Nmn, vswf1ds.Mmn), ndS)
									+ rap
											* dotProduct(
													crossProduct(vswf2ds.Mmn,
															vswf1ds.Nmn), ndS));
					c++;
				}
				l++;
			}
		}
	}
	for (int i = 0; i < 2 * N + 1; i++) {
		gaussj(Q312[i], Q311[i]);
		Matrix X = Q311[i];
		T[i] = Q11[i] * X;
	}
	delete[] Q11;
	delete[] Q311;
	delete[] Q312;
	return;
}

Matrix AxisSymParticle::getT() {
	int NM = (N + 1) * (N + 1) - 1;
	Matrix TT(2 * NM, 2 * NM, 0.);
	int nbr, l, c;
	for (int m = -N; m <= N; m++) {
		l = 0;
		nbr = N - abs(m) + 1;
		if (abs(m) == 0)
			nbr--;
		for (int nl = abs(m); nl <= N; nl++) {
			if (nl == 0)
				nl++;
			c = 0;
			for (int nc = abs(m); nc <= N; nc++) {
				if (nc == 0)
					nc++;
				TT.t[ii(nl, m)][ii(nc, m)] = T[m + N].t[l][c];
				TT.t[ii(nl, m)][ii(nc, m) + NM] = T[m + N].t[l][c + nbr];
				TT.t[ii(nl, m) + NM][ii(nc, m)] = T[m + N].t[l + nbr][c];
				TT.t[ii(nl, m) + NM][ii(nc, m) + NM] = T[m + N].t[l + nbr][c
						+ nbr];
				c++;
			}
			l++;
		}
	}
	return TT;
}

Matrix AxisSymParticle::getR() {
	int NM = (N + 1) * (N + 1) - 1;
	Matrix RR(2 * NM, 2 * NM, 0.);
	int nbr, l, c;
	for (int m = -N; m <= N; m++) {
		l = 0;
		nbr = N - abs(m) + 1;
		if (abs(m) == 0)
			nbr--;
		for (int nl = abs(m); nl <= N; nl++) {
			if (nl == 0)
				nl++;
			c = 0;
			for (int nc = abs(m); nc <= N; nc++) {
				if (nc == 0)
					nc++;
				RR.t[ii(nl, m)][ii(nc, m)] = R[m + N].t[l][c];
				RR.t[ii(nl, m)][ii(nc, m) + NM] = R[m + N].t[l][c + nbr];
				RR.t[ii(nl, m) + NM][ii(nc, m)] = R[m + N].t[l + nbr][c];
				RR.t[ii(nl, m) + NM][ii(nc, m) + NM] = R[m + N].t[l + nbr][c
						+ nbr];
				c++;
			}
			l++;
		}
	}
	return RR;
}

Vector AxisSymParticle::multT(const Vector &v) {
	int NM = (N + 1) * (N + 1) - 1;
	Vector vv(2 * NM, 0.);
	int nbr, l, c;
	for (int m = -N; m <= N; m++) {
		l = 0;
		nbr = N - abs(m) + 1;
		if (abs(m) == 0)
			nbr--;
		for (int nl = abs(m); nl <= N; nl++) {
			if (nl == 0)
				nl++;
			c = 0;
			for (int nc = abs(m); nc <= N; nc++) {
				if (nc == 0)
					nc++;
				vv.t[ii(nl, m)] += T[m + N].t[l][c] * v.t[ii(nc, m)]
						+ T[m + N].t[l][c + nbr] * v.t[ii(nc, m) + NM];
				vv.t[ii(nl, m) + NM] += T[m + N].t[l + nbr][c] * v.t[ii(nc, m)]
						+ T[m + N].t[l + nbr][c + nbr] * v.t[ii(nc, m) + NM];
				c++;
			}
			l++;
		}
	}
	return vv;
}
