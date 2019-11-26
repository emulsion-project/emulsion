#include "MultipleScattering.h"
#include <iostream>
using namespace std;

MultipleScattering::MultipleScattering(const string &fname) {
	string str;
	ifstream ifs(fname.c_str());
	int mat, typ, inc;
	RVector geom;
	ifs >> str >> N >> str >> Ng >> str >> methodTMatrix >> str >> solver >> str
			>> inc >> str >> lambda;
	ifs >> str >> xc >> yc >> zc >> str >> Np;
	ifs >> str >> periodic >> nbrDim >> str >> Npx >> Npy >> Npz >> str >> dx
			>> dy >> dz;
	NM = (N + 1) * (N + 1) - 1;

	particles = new AxisSymParticle[Np];
	ml.setParameters("Multilayer/multilayerConfig.txt", lambda);

	x = new double[Np];
	y = new double[Np];
	z = new double[Np];

	positions = RMatrix(Np, 3);
	ifs >> str;
	for (int i = 0; i < Np; i++) {
		ifs >> mat >> typ;
		ifs >> x[i] >> y[i] >> z[i];
		positions.t[i][0] = x[i];
		positions.t[i][1] = y[i];
		positions.t[i][2] = z[i];
		if (typ == 0) {
			geom = RVector(1);
			ifs >> geom.t[0];
		} else if (typ == 1 || typ == 2) {
			geom = RVector(2);
			ifs >> geom.t[0] >> geom.t[1];
		}
		particles[i].setParticleProperties(typ, mat,
				ml.mediums[ml.getMedium(z[i])], lambda, geom);
		particles[i].setTMatrixParameters(N, Ng, methodTMatrix);
	}

	lay = ml.getMedium(z[0]);
	if (ml.nbrMediums == 1)
		z0 = 0.;
	else if (lay == 0)
		z0 = ml.Zi.t[0];
	else if (lay == ml.nbrMediums - 1)
		z0 = ml.Zi.t[ml.nbrMediums - 2];
	else
		z0 = (ml.Zi.t[lay - 1] + ml.Zi.t[lay]) / 2.;
	if (ml.nbrMediums == 1)
		z0 = 0.;
	incidence.setParameters(inc, lambda);

	//refl.setConstants(lambda,N,Ng);
	refl1.setConstants(lambda, N, Ng);

	matricesFFT1D = NULL;
}

void MultipleScattering::setN(const int &n) {
	N = n;
	NM = (N + 1) * (N + 1) - 1;
	//refl.setConstants(lambda,N,Ng);
	refl1.setConstants(lambda, N, Ng);
	for (int i = 0; i < Np; i++)
		particles[i].setTMatrixParameters(N, Ng, methodTMatrix);
	return;
}

void MultipleScattering::setLambda(const double &l) {
	lambda = l;
	incidence.setParameters(incidence.incidence, lambda);
	ml.setLambda(lambda);
	//if (ml.nbrMediums>1) refl.setConstants(lambda,N,Ng);
	if (ml.nbrMediums > 1)
		refl1.setConstants(lambda, N, Ng);
	for (int i = 0; i < Np; i++)
		particles[i].setLambda(lambda);
	for (int i = 0; i < Np; i++)
		particles[i].setParticleProperties(particles[i].type,
				particles[i].material, ml.mediums[ml.getMedium(z[i])], lambda,
				particles[i].geom);
	for (int i = 0; i < Np; i++)
		particles[i].setTMatrixParameters(N, Ng, methodTMatrix);
	return;
}

void MultipleScattering::setElectronEnergy(const double &ee) {
	setLambda(2. * piDbl * c * hbar / (ee * e));
	return;
}

void MultipleScattering::calcAmni() {
	Vector3d pos;
	Amni = Vector(2 * NM * Np);
	for (int n = 0; n < Np; n++) {
		pos.t[0] = x[n];
		pos.t[1] = y[n];
		pos.t[2] = z[n];
		incidence.setPlaneWaveCoeffs(N, pos);
		for (int i = 0; i < 2 * NM; i++)
			Amni.t[i + 2 * n * NM] = incidence.Amn.t[i];
	}
	return;
}

void MultipleScattering::calcTAmni() {
	Vector3d pos;
	Amni = Vector(2 * NM * Np, 0.);
	Vector A(2 * NM);
	for (int n = 0; n < Np; n++) {
		pos.t[0] = x[n];
		pos.t[1] = y[n];
		pos.t[2] = z[n];
		incidence.setPlaneWaveCoeffs(N, pos);
		A = particles[n].multT(incidence.Amn);
		for (int i = 0; i < 2 * NM; i++)
			Amni.t[i + 2 * n * NM] = A.t[i];
	}
	return;
}

void MultipleScattering::calcAmniElectron() {
	Amni = Vector(2 * NM * Np);
	for (int n = 0; n < Np; n++) {
		incidence.setElectron(N, incidence.E0, incidence.E, x[n], y[n]);
		for (int i = 0; i < 2 * NM; i++)
			Amni.t[i + 2 * n * NM] = incidence.Amn.t[i];
	}
	return;
}

void MultipleScattering::calcTAmniElectron() {
	Amni = Vector(2 * NM * Np);
	Vector A(2 * NM);
	for (int n = 0; n < Np; n++) {
		incidence.setElectron(N, incidence.E0, incidence.E, x[n], y[n]);
		for (int i = 0; i < 2 * NM; i++)
			A.t[i] = incidence.Amn.t[i];
		A = particles[n].multT(A);
		for (int i = 0; i < 2 * NM; i++)
			Amni.t[i + 2 * n * NM] = A.t[i];
	}
	return;
}

void MultipleScattering::calcPmn() {
	Pmn = Vector(2 * NM, 0.);
	Vector A(2 * NM);
	for (int i = 0; i < Np; i++) {
		for (int ii = 0; ii < 2 * NM; ii++)
			A.t[ii] = Pmni.t[ii + 2 * i * NM];
		trans.setParameters(1, false, N, particles[0].k1, -x[i], -y[i], -z[i]);
		A = trans.translateT(A, 0);
		for (int ii = 0; ii < 2 * NM; ii++)
			Pmn.t[ii] += A.t[ii];
	}
}

Vector MultipleScattering::calcPmni0() {
	Vector Pmni0(Np * 2 * NM);
	Vector A(2 * NM);
	for (int i = 0; i < Np; i++) {
		for (int ii = 0; ii < 2 * NM; ii++)
			A.t[ii] = Pmni.t[ii + 2 * i * NM];
		trans.setParameters(1, false, N, particles[0].k1, -x[i], -y[i], -z[i]);
		A = trans.translateT(A, 0);
		for (int ii = 0; ii < 2 * NM; ii++)
			Pmni0.t[ii + 2 * i * NM] = A.t[ii];
	}
	return Pmni0;
}

void MultipleScattering::solve() {
	if (incidence.type == 0)
		calcTAmni();
	else if (incidence.type == 1)
		calcTAmniElectron();
	if ((Np == 1) && (ml.nbrMediums == 1))
		Pmni = Amni;
	else if (solver == 0)
		Pmni = calcDirect() * Amni;
	else if (solver == 1 || solver == 2) {
		if (periodic == 1) {
			if (nbrDim == 1)
				constructMatrix1D();
			else if (nbrDim == 2)
				constructMatrix2D();
		}
		if (ml.nbrMediums > 1)
			refl1.setPositions(positions);
		gmres(Amni, 1e-5);
		Pmni = Amni;

		if (periodic == 1) {
			delete[] matricesFFT1D;
		}
	}
	calcAmni();
	Vector cmn(2 * NM);
	Cmni = Vector(2 * NM * Np, 0.);
	for (int i = 0; i < Np; i++) {
		for (int j = 0; j < 2 * NM; j++)
			cmn.t[j] = Amni.t[j + i * 2 * NM];
		cmn = particles[i].getR() * cmn;
		for (int j = 0; j < 2 * NM; j++)
			Cmni.t[j + i * 2 * NM] = cmn.t[j];
	}
}

Matrix MultipleScattering::calcDirect() {
	Matrix M(2 * NM * Np, 2 * NM * Np, 0.);
	Matrix T;
	for (int i = 0; i < Np; i++) {
		for (int j = i; j < Np; j++) {
			if (i == j)
				for (int ii = 0; ii < 2 * NM; ii++)
					M.t[ii + i * 2 * NM][ii + i * 2 * NM] = 1.;
			else {
				trans.setParameters(3, true, N, particles[0].k1, x[i] - x[j],
						y[i] - y[j], z[i] - z[j]);
				T = particles[i].getT() * trans.getMatrixT(0);
				for (int ii = 0; ii < 2 * NM; ii++)
					for (int jj = 0; jj < 2 * NM; jj++)
						M.t[ii + i * 2 * NM][jj + j * 2 * NM] = -T.t[ii][jj];
				T = particles[j].getT() * trans.getMatrixT(1);
				for (int ii = 0; ii < 2 * NM; ii++)
					for (int jj = 0; jj < 2 * NM; jj++)
						M.t[ii + j * 2 * NM][jj + i * 2 * NM] = -T.t[ii][jj];
			}
		}
	}
	invGaussj(M);
	return M;
}

Vector MultipleScattering::calcDirectGMRES(const Vector &B) {
	if ((Np == 1) && (ml.nbrMediums == 1))
		return B;
	Vector A(2 * NM), X(2 * NM * Np, 0.);
	int ii;
	for (int i = 0; i < Np; i++) {
		for (int j = i; j < Np; j++) {
			if (i == j) {
				if (ml.nbrMediums > 1) {
					//refl.setParameters(0.,0.,z[i],z[i],true);
					refl1.calcRMatrices(i, i);
					for (ii = 0; ii < 2 * NM; ii++)
						A.t[ii] = B.t[ii + 2 * NM * j];
					A = transpose(refl1.R) * A;
					for (ii = 0; ii < 2 * NM; ii++)
						X.t[ii + 2 * NM * i] -= A.t[ii];
				}
			} else {
				if (ml.nbrMediums == 1) {
					trans.setParameters(3, false, N, particles[0].k1,
							x[i] - x[j], y[i] - y[j], z[i] - z[j]);
					for (ii = 0; ii < 2 * NM; ii++)
						A.t[ii] = B.t[ii + 2 * NM * j];
					A = trans.translateT(A, 0);
					for (ii = 0; ii < 2 * NM; ii++) {
						X.t[ii + 2 * NM * i] -= A.t[ii];
						A.t[ii] = B.t[ii + 2 * NM * i];
					}
					A = trans.translateT(A, 1);
					for (ii = 0; ii < 2 * NM; ii++)
						X.t[ii + 2 * NM * j] -= A.t[ii];
				} else {
					if (ml.getMedium(z[i]) == ml.getMedium(z[j])) {
						trans.setParameters(3, false, N, particles[0].k1,
								x[i] - x[j], y[i] - y[j], z[i] - z[j]);
						for (ii = 0; ii < 2 * NM; ii++)
							A.t[ii] = B.t[ii + 2 * NM * j];
						A = trans.translateT(A, 0);
						for (ii = 0; ii < 2 * NM; ii++) {
							X.t[ii + 2 * NM * i] -= A.t[ii];
							A.t[ii] = B.t[ii + 2 * NM * i];
						}
						A = trans.translateT(A, 1);
						for (ii = 0; ii < 2 * NM; ii++)
							X.t[ii + 2 * NM * j] -= A.t[ii];
					}
					//refl.setParameters(x[i]-x[j],y[i]-y[j],z[j],z[i],false);
					refl1.calcRMatrices(i, j);
					for (ii = 0; ii < 2 * NM; ii++)
						A.t[ii] = B.t[ii + 2 * NM * j];
					A = transpose(refl1.R) * A;
					for (ii = 0; ii < 2 * NM; ii++)
						X.t[ii + 2 * NM * i] -= A.t[ii];
					//refl.setParameters(x[j]-x[i],y[j]-y[i],z[i],z[j],false);
					for (ii = 0; ii < 2 * NM; ii++)
						A.t[ii] = B.t[ii + 2 * NM * i];
					A = transpose(refl1.Ri) * A;
					for (ii = 0; ii < 2 * NM; ii++)
						X.t[ii + 2 * NM * j] -= A.t[ii];
				}
			}
		}
	}
	for (int i = 0; i < Np; i++) {
		for (ii = 0; ii < 2 * NM; ii++)
			A.t[ii] = X.t[ii + 2 * NM * i];
		A = particles[i].multT(A);
		for (ii = 0; ii < 2 * NM; ii++)
			X.t[ii + 2 * NM * i] = A.t[ii];
	}
	for (ii = 0; ii < 2 * NM * Np; ii++)
		X.t[ii] += B.t[ii];
	return X;
}

Vector MultipleScattering::calcFFT1D(const Vector &B) {
	Vector X(2 * NM * Npx);
	Vector A(2 * NM);
	Vector AA(nMatFFT);
	Vector *X1;
	X1 = new Vector[nMatFFT];
	for (int i = 0; i < nMatFFT; i++)
		X1[i] = Vector(2 * NM, 0.);
	for (int i = 0; i < Npx; i++)
		for (int j = 0; j < 2 * NM; j++)
			X1[i].t[j] = B.t[j + 2 * NM * i];
	for (int i = 0; i < 2 * NM; i++) {
		for (int j = 0; j < nMatFFT; j++)
			AA.t[j] = X1[j].t[i];
		AA = fft(AA, 1);
		for (int j = 0; j < nMatFFT; j++)
			X1[j].t[i] = AA.t[j];
	}
	for (int i = 0; i < nMatFFT; i++)
		X1[i] = matricesFFT1D[i] * X1[i];
	for (int i = 0; i < 2 * NM; i++) {
		for (int j = 0; j < nMatFFT; j++)
			AA.t[j] = X1[j].t[i];
		AA = fft(AA, -1);
		for (int j = 0; j < Npx; j++)
			X.t[i + j * 2 * NM] = AA.t[j];
	}
	for (int i = 0; i < Npx; i++) {
		for (int ii = 0; ii < 2 * NM; ii++)
			A.t[ii] = X.t[ii + 2 * NM * i];
		A = particles[i].multT(A);
		for (int ii = 0; ii < 2 * NM; ii++)
			X.t[ii + 2 * NM * i] = A.t[ii];
	}
	for (int ii = 0; ii < 2 * NM * Np; ii++)
		X.t[ii] = B.t[ii] - X.t[ii];
	return X;
}

void MultipleScattering::constructMatrix1D() {
	nMatFFT = 2;
	for (int i = 1; 2 * Npx > nMatFFT; i++)
		nMatFFT += nMatFFT;
	matricesFFT1D = new Matrix[nMatFFT];
	for (int i = 0; i < nMatFFT; i++)
		matricesFFT1D[i] = Matrix(2 * NM, 2 * NM);
	Matrix *matrices;
	matrices = new Matrix[nMatFFT];
	for (int i = 0; i < nMatFFT; i++)
		matrices[i] = Matrix(2 * NM, 2 * NM, 0.);
	for (int i = 1; i < Npx; i++) {
		trans.setParameters(3, true, N, particles[0].k1, x[0] - x[i], 0., 0.);
		matrices[i] = trans.getMatrixT(0);
		matrices[nMatFFT - i] = trans.getMatrixT(1);
	}
	Vector v(nMatFFT);
	for (int i = 0; i < 2 * NM; i++) {
		for (int j = 0; j < 2 * NM; j++) {
			for (int k = 0; k < nMatFFT; k++)
				v.t[k] = matrices[k].t[i][j];
			v = fft(v, 1);
			for (int k = 0; k < nMatFFT; k++)
				matricesFFT1D[k].t[i][j] = v.t[k];
		}
	}
	delete[] matrices;
	return;
}

Vector MultipleScattering::calcFFT2D(const Vector &B) {
	Vector X(2 * NM * Npx * Npy);
	Vector A(2 * NM);
	Vector AAx(nMatFFTx), AAy(nMatFFTy);
	Vector *X1;
	X1 = new Vector[nMatFFTx * nMatFFTy];
	for (int i = 0; i < nMatFFTx * nMatFFTy; i++)
		X1[i] = Vector(2 * NM, 0.);
	for (int i = 0; i < Npx; i++)
		for (int ii = 0; ii < Npy; ii++)
			for (int j = 0; j < 2 * NM; j++)
				X1[i * nMatFFTy + ii].t[j] = B.t[j + 2 * NM * (i * Npy + ii)];
	for (int i = 0; i < 2 * NM; i++) {
		for (int ii = 0; ii < Npy; ii++) {
			for (int j = 0; j < nMatFFTx; j++)
				AAx.t[j] = X1[j * nMatFFTy + ii].t[i];
			AAx = fft(AAx, 1);
			for (int j = 0; j < nMatFFTx; j++)
				X1[j * nMatFFTy + ii].t[i] = AAx.t[j];
		}
		for (int ii = 0; ii < nMatFFTx; ii++) {
			for (int j = 0; j < nMatFFTy; j++)
				AAy.t[j] = X1[ii * nMatFFTy + j].t[i];
			AAy = fft(AAy, 1);
			for (int j = 0; j < nMatFFTy; j++)
				X1[ii * nMatFFTy + j].t[i] = AAy.t[j];
		}
	}
	for (int i = 0; i < nMatFFTx * nMatFFTy; i++)
		X1[i] = matricesFFT1D[i] * X1[i];
	for (int i = 0; i < 2 * NM; i++) {
		for (int ii = 0; ii < nMatFFTx; ii++) {
			for (int j = 0; j < nMatFFTy; j++)
				AAy.t[j] = X1[ii * nMatFFTy + j].t[i];
			AAy = fft(AAy, -1);
			for (int j = 0; j < Npy; j++)
				X1[ii * nMatFFTy + j].t[i] = AAy.t[j];
		}
		for (int ii = 0; ii < Npy; ii++) {
			for (int j = 0; j < nMatFFTx; j++)
				AAx.t[j] = X1[j * nMatFFTy + ii].t[i];
			AAx = fft(AAx, -1);
			for (int j = 0; j < Npx; j++)
				X.t[i + (j * Npy + ii) * 2 * NM] = AAx.t[j];
		}
	}
	for (int i = 0; i < Npx * Npy; i++) {
		for (int ii = 0; ii < 2 * NM; ii++)
			A.t[ii] = X.t[ii + 2 * NM * i];
		A = particles[i].multT(A);
		for (int ii = 0; ii < 2 * NM; ii++)
			X.t[ii + 2 * NM * i] = A.t[ii];
	}
	for (int ii = 0; ii < 2 * NM * Npx * Npy; ii++)
		X.t[ii] = B.t[ii] - X.t[ii];
	delete[] X1;
	return X;
}

void MultipleScattering::constructMatrix2D() {
	//Npx = Npy = 4;
	nMatFFTx = nMatFFTy = 2;
	for (int i = 1; 2 * Npx > nMatFFTx; i++)
		nMatFFTx += nMatFFTx;
	for (int i = 1; 2 * Npy > nMatFFTy; i++)
		nMatFFTy += nMatFFTy;
	matricesFFT1D = new Matrix[nMatFFTx * nMatFFTy];
	for (int i = 0; i < nMatFFTx * nMatFFTy; i++)
		matricesFFT1D[i] = Matrix(2 * NM, 2 * NM);
	Matrix *matrices;
	matrices = new Matrix[nMatFFTx * nMatFFTy];
	for (int i = 0; i < nMatFFTx * nMatFFTy; i++)
		matrices[i] = Matrix(2 * NM, 2 * NM, 0.);
	for (int i = 1; i < Npx; i++) {
		trans.setParameters(3, true, N, particles[0].k1, x[0] - x[i * Npy], 0.,
				0.);
		matrices[i * nMatFFTy] = trans.getMatrixT(0);
		matrices[(nMatFFTx - i) * nMatFFTy] = trans.getMatrixT(1);
	}
	for (int i = 1; i < Npy; i++) {
		trans.setParameters(3, true, N, particles[0].k1, 0., y[0] - y[i], 0.);
		matrices[i] = trans.getMatrixT(0);
		matrices[nMatFFTy - i] = trans.getMatrixT(1);
	}
	for (int i = 1; i < Npx; i++) {
		for (int j = 1; j < Npy; j++) {
			trans.setParameters(3, true, N, particles[0].k1, x[0] - x[i * Npy],
					y[0] - y[j], 0.);
			matrices[i * nMatFFTy + j] = trans.getMatrixT(0);
			matrices[(nMatFFTx - i + 1) * nMatFFTy - j] = trans.getMatrixT(1);
			trans.setParameters(3, true, N, particles[0].k1, x[0] - x[i * Npy],
					y[j] - y[0], 0.);
			matrices[(i + 1) * nMatFFTy - j] = trans.getMatrixT(0);
			matrices[(nMatFFTx - i) * nMatFFTy + j] = trans.getMatrixT(1);
		}
	}
	Vector vx(nMatFFTx), vy(nMatFFTy);
	for (int i = 0; i < 2 * NM; i++) {
		for (int j = 0; j < 2 * NM; j++) {
			for (int ii = 0; ii < Npy; ii++) {
				for (int jj = 0; jj < nMatFFTx; jj++)
					vx.t[jj] = matrices[jj * nMatFFTy + ii].t[i][j];
				vx = fft(vx, 1);
				for (int jj = 0; jj < nMatFFTx; jj++)
					matricesFFT1D[jj * nMatFFTy + ii].t[i][j] = vx.t[jj];
			}
			for (int ii = nMatFFTy - Npy + 1; ii < nMatFFTy; ii++) {
				for (int jj = 0; jj < nMatFFTx; jj++)
					vx.t[jj] = matrices[jj * nMatFFTy + ii].t[i][j];
				vx = fft(vx, 1);
				for (int jj = 0; jj < nMatFFTx; jj++)
					matricesFFT1D[jj * nMatFFTy + ii].t[i][j] = vx.t[jj];
			}
			for (int ii = 0; ii < nMatFFTx; ii++) {
				for (int jj = 0; jj < nMatFFTy; jj++)
					vy.t[jj] = matricesFFT1D[ii * nMatFFTy + jj].t[i][j];
				vy = fft(vy, 1);
				for (int jj = 0; jj < nMatFFTy; jj++)
					matricesFFT1D[ii * nMatFFTy + jj].t[i][j] = vy.t[jj];
			}
		}
	}
	delete[] matrices;
	return;
}

void MultipleScattering::gmres(Vector &X, const double &tol) {
	const int kmax = 1000;
	int i, k, it, itt = 0, n = X.n;
	double res;
	Complex *h1, *h2, aux, ss[kmax], cc[kmax], gg[kmax + 1], y[kmax];
	Vector XX, B, *V[kmax + 1], *H[kmax + 1];
	XX = B = X;
	Matrix M;
	do {
		it = 0;
		for (i = 0; i < n; i++)
			X.t[i] = 0.;
		for (i = 0; i < kmax + 1; i++)
			gg[i] = 0.;
		V[0] = new Vector(n);
		gg[0] = res = norm(XX);
		if (res > tol)
			res = 1. / res;
		for (i = 0; i < n; i++)
			V[0]->t[i] = res * XX.t[i];
		do {
			V[it + 1] = new Vector(n);
			H[it] = new Vector(it + 2);
			if (solver == 1)
				*V[it + 1] = calcDirectGMRES(*V[it]);
			else if (solver == 2) {
				if (nbrDim == 1)
					*V[it + 1] = calcFFT1D(*V[it]);
				else if (nbrDim == 2)
					*V[it + 1] = calcFFT2D(*V[it]);
			}
			for (k = 0; k < it + 1; k++) {
				H[it]->t[k] = aux = (*V[k]) * (*V[it + 1]);
				for (i = 0; i < n; i++)
					V[it + 1]->t[i] -= aux * V[k]->t[i];
			}
			H[it]->t[it + 1] = res = norm(*V[it + 1]);
			if (fabs(res) > 1e-10) {
				res = 1. / res;
				for (i = 0; i < n; i++)
					V[it + 1]->t[i] *= res;
			}
			h1 = h2 = H[it]->t;
			h2++;
			for (k = 0; k < it; k++, h1++, h2++) {
				*h1 = cc[k] * (aux = *h1) + ss[k] * (*h2);
				*h2 = -conj(ss[k]) * aux + cc[k] * (*h2);
			}
			h1 = h2 = H[it]->t + it;
			h2++;
			if (abs(*h2) < 1e-10) {
				ss[it] = 0.;
				cc[it] = 1.;
			} else if (abs(*h1) < 1e-10) {
				ss[it] = conj(*h2) / abs(*h2);
				cc[it] = 0.;
			} else {
				res = 1. / sqrt(norm(*h1) + norm(*h2));
				cc[it] = res * abs(*h1);
				ss[it] = res * conj(*h2) * (*h1) / abs(*h1);
			}
			*h1 = cc[it] * (aux = *h1) + ss[it] * (*h2);
			*h2 = -conj(ss[it]) * aux + cc[it] * (*h2);
			gg[it + 1] = -gg[it] * conj(ss[it]);
			gg[it] *= cc[it];
			res = abs(gg[it + 1]);
			//cout << it << " residual " << res << endl;
			it++;
			itt++;
		} while ((it < kmax) && (res > tol));
		for (k = it - 1; k >= 0; k--) {
			y[k] = gg[k];
			for (i = it - 1; i > k; i--)
				y[k] -= y[i] * H[i]->t[k];
			y[k] /= H[k]->t[k];
		}
		for (k = 0; k < it; k++)
			for (i = 0; i < n; i++)
				X.t[i] += y[k] * V[k]->t[i];
		for (k = 0; k < it; k++) {
			delete V[k];
			delete H[k];
		}
		delete V[it];
		if (solver == 1)
			XX = B - calcDirectGMRES(X);
		else if (solver == 2) {
			if (nbrDim == 1)
				XX = B - calcFFT1D(X);
			else if (nbrDim == 2)
				XX = B - calcFFT2D(X);
		}
	} while (res > tol);
}

double MultipleScattering::calcCsca(const bool &normalize, const double &nor) {
	calcPmn();
	double Csca = 0.;
	double mult = piDbl / (real(particles[0].k1 * particles[0].k1));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < 2 * NM; i++)
		Csca += norm(Pmn.t[i]);
	return Csca * mult;
}

RVector MultipleScattering::calcCext(const bool &normalize, const double &nor) {
	calcAmni();
	RVector Cext(Np, 0.);
	double mult = piDbl / (real(particles[0].k1 * particles[0].k1));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++)
		for (int j = 0; j < 2 * NM; j++)
			Cext.t[i] -= real(
					Pmni.t[j + 2 * NM * i] * conj(Amni.t[j + 2 * NM * i]))
					* mult;
	return Cext;
}

RVector MultipleScattering::calcCextOT(const bool &normalize,
		const double &nor) {
	RVector Cext(8, 0.);
	RVector inc(4);
	Vector3d Ep, Em;
	Vector is, ip;
	is = incidence.aouts;
	ip = incidence.aoutp;
	Ep = calcEscaInf(
			acos(
					incidence.ml.calcKz(ml.nbrMediums - 1, true)
							/ ml.ki.t[ml.nbrMediums - 1]).re, incidence.alpha0);
	//cout << "ext\t" << acos(incidence.ml.calcKz(ml.nbrMediums-1,true)/ml.ki.t[ml.nbrMediums-1]).re*180./piDbl << "\t" << incidence.alpha0 << "\t" << Ep << endl;
	Em = calcEscaInf(piDbl - acos(incidence.ml.calcKz(0, true) / ml.ki.t[0]).re,
			incidence.alpha0);
	double multp = 4. * piDbl
			/ (ml.ki.t[ml.nbrMediums - 1].re * (norm(ip.t[1]) + norm(is.t[1])));
	double multm = 4. * piDbl
			/ (ml.ki.t[0].re * (norm(ip.t[0]) + norm(is.t[0])));
	if (normalize)
		multp /= piDbl * nor * nor;
	if (normalize)
		multm /= piDbl * nor * nor;
	Cext.t[0] = imag(
			multp * (Ep.t[1] * conj(ip.t[1]) + Ep.t[2] * conj(is.t[1])));
	Cext.t[1] = imag(
			multm * (Em.t[1] * conj(ip.t[0]) + Em.t[2] * conj(is.t[0])));
	Cext.t[2] = phase(Ep.t[2]);
	Cext.t[3] = phase(Em.t[2]);
	Cext.t[4] = phase(is.t[1]);
	Cext.t[5] = phase(is.t[0]);
	Cext.t[6] = dotProduct(Ep, Ep).re;
	Cext.t[7] = dotProduct(Em, Em).re;

	return Cext;
}

RVector MultipleScattering::calcScaML(const bool &normalize,
		const double &nor) {
	RVector Csca(2, 0.);
	Vector3d Ep, Em;
	Vector is, ip;
	RMatrix gaussAlpha = gaussLegQuad(Ng / 2, 0., 2. * piDbl);
	RMatrix gaussBeta1 = gaussLegQuad(Ng / 2, 0., piDbl / 2.);
	RMatrix gaussBeta2 = gaussLegQuad(Ng / 2, piDbl / 2., piDbl);
	for (int ka = 0; ka < Ng / 4; ka++) {
		for (int kb = 0; kb < Ng / 4; kb++) {
			Ep = calcEscaInf(gaussBeta1.t[kb][0], gaussAlpha.t[ka][0]);
			Em = calcEscaInf(gaussBeta2.t[kb][0], gaussAlpha.t[ka][0]);
			Csca.t[0] += gaussAlpha.t[ka][1] * gaussBeta1.t[kb][1]
					* dotProduct(Ep, Ep).re;
			Csca.t[1] += gaussAlpha.t[ka][1] * gaussBeta1.t[kb][1]
					* dotProduct(Em, Em).re;
		}
	}
	return Csca;
}

double MultipleScattering::calcCscaN(const bool &normalize, const double &nor,
		const int &nn, const int &pol) {
	calcPmn();
	double Csca = 0.;
	double mult = piDbl / (real(particles[0].k1 * particles[0].k1));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = -nn; i <= nn; i++)
		Csca += norm(Pmn.t[ii(nn, i) + pol * NM]);
	return Csca * mult;
}

RVector MultipleScattering::calcCextN(const bool &normalize, const double &nor,
		const int &nn, const int &pol) {
	calcAmni();
	Vector Pmni0 = calcPmni0();
	RVector Cext(Np, 0.);
	double mult = piDbl / (real(particles[0].k1 * particles[0].k1));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++)
		for (int j = -nn; j <= nn; j++)
			Cext.t[i] -= real(
					Pmni0.t[ii(nn, j) + pol * NM + 2 * NM * i]
							* conj(Amni.t[ii(nn, j) + pol * NM])) * mult;
	return Cext;
}

RVector MultipleScattering::calcCabsN(const bool &normalize, const double &nor,
		const int &nn, const int &pol) {
	calcAmni();
	RVector Cabs(Np, 0.);
	Vector A(2 * NM);
	double mult = piDbl
			/ (real(particles[0].k1 * particles[0].k1)
					* (norm(incidence.Ea) + norm(incidence.Eb)));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++) {
		A = Vector(2 * NM, 1.);
		A = particles[i].getT() * A;
		for (int j = -nn; j <= nn; j++)
			Cabs.t[i] -= mult
					* (real(1. / A.t[ii(nn, j) + pol * NM] + 1.)
							* norm(Pmni.t[ii(nn, j) + pol * NM + 2 * NM * i]));
	}
	return Cabs;
}

RVector MultipleScattering::calcCabs(const bool &normalize, const double &nor) {
	calcAmni();
	RVector Cabs(Np, 0.);
	Vector A(2 * NM);
	double mult = piDbl
			/ (real(particles[0].k1 * particles[0].k1)
					* (norm(incidence.Ea) + norm(incidence.Eb)));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++) {
		A = Vector(2 * NM, 1.);
		A = particles[i].getT() * A;
		for (int j = 0; j < 2 * NM; j++)
			Cabs.t[i] -= mult
					* (real(1. / A.t[j] + 1.) * norm(Pmni.t[j + 2 * NM * i]));
	}
	return Cabs;
}

Vector MultipleScattering::calcComplexCext(const bool &normalize,
		const double &nor) {
	calcAmni();
	Vector Cext(Np, 0.);
	double mult = piDbl
			/ (real(particles[0].k1 * particles[0].k1)
					* (norm(incidence.Ea) + norm(incidence.Eb)));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++)
		for (int j = 0; j < 2 * NM; j++)
			Cext.t[i] -= mult * conj(incidence.Amn.t[j])
					* Pmni.t[j + 2 * NM * i] * i_;
	return Cext;
}

Vector MultipleScattering::calcComplexCextN(const bool &normalize,
		const double &nor, const int &nn, const int &pol) {
	calcAmni();
	//Vector Pmni0 = calcPmni0();
	Vector Cext(Np, 0.);
	double mult = piDbl / (real(particles[0].k1 * particles[0].k1));
	if (normalize)
		mult /= piDbl * nor * nor;
	for (int i = 0; i < Np; i++)
		for (int j = -nn; j <= nn; j++)
			Cext.t[i] -= (Pmni.t[ii(nn, j) + pol * NM + 2 * NM * i]
					* conj(incidence.Amn.t[ii(nn, j) + pol * NM])) * mult;
	return Cext;
}

double MultipleScattering::calcEEL() {
	calcAmniElectron();
	double xx = 0.;
	double mult;
	int ii = 0;
	double omega = incidence.omega;
	for (int i = 0; i < Np; i++) {
		ii = 0;
		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				mult = c * c * double(n * n + n) / (4. * piDbl * piDbl * omega);
				xx -= mult
						* real(
								conj(Amni.t[ii + 2 * NM * i])
										* Pmni.t[ii + 2 * NM * i]);
				xx -= mult
						* real(
								conj(Amni.t[ii + 2 * NM * i + NM])
										* Pmni.t[ii + 2 * NM * i + NM]);
				ii++;
			}
		}
	}
	return xx;
}

double MultipleScattering::calcCL() {
	//calcAmniElectron();
	calcPmn();
	double xx = 0.;
	double mult;
	//int ii=0;
	for (int n = 1; n <= N; n++) {
		for (int m = -n; m <= n; m++) {
			mult = double(n * n + n)
					/ (4. * piDbl * piDbl * particles[0].k1.re
							* particles[0].k1.re * particles[0].k1.re);
			xx += mult * (norm(Pmn.t[ii(n, m)]));
			xx += mult * (norm(Pmn.t[ii(n, m) + NM]));
			//xx += mult*double(m)*v0*real(conj(Amnplus(m,n,vc))*powi(n)*Pmni.t[ii+2*NM*i]*exp(i_*double(m)*phi0+i_*omega*z0/v0));
			//xx += mult*c/(2.*gamma)*real(conj(Bmnplus(m,n,vc))*powi(n)*Pmni.t[ii+2*NM*i+NM]*exp(i_*double(m)*phi0+i_*omega*z0/v0));
			//cout << n << "\t" << m << "\t" << xx << "\t" << imag(conj(Amnplus(m,n,vc))*powi(n)*Pmni.t[ii+2*NM*i]*exp(i_*double(m)*phi0+i_*omega*z0/v0)) << "\t" << imag(conj(Bmnplus(m,n,vc))*powi(n)*Pmni.t[ii+2*NM*i+NM]*exp(i_*double(m)*phi0+i_*omega*z0/v0)) << endl;
			//ii++;
		}
	}
	return xx;
}

double MultipleScattering::calcEEL(const int &n, const int &m,
		const int &mode) {
	calcAmniElectron();
	double xx = 0.;
	double mult;
	double omega = incidence.omega;
	for (int i = 0; i < Np; i++) {
		mult = c * c * double(n * n + n) / (4. * piDbl * piDbl * omega);
		xx -= mult
				* real(
						conj(Amni.t[ii(n, m) + 2 * NM * i + mode * NM])
								* Pmni.t[ii(n, m) + 2 * NM * i + mode * NM]);
	}
	return xx;
}

Vector3d MultipleScattering::calcEsca(const double &xx, const double &yy,
		const double &zz) {
	Vector3d E(0.);
	Vector3d Mmn, Nmn;
	double X, Y, Z, R, T, P;
	int ii;
	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		if (R < particles[i].calcR(Z / R) * 0.999)
			return calcEint(i, xx, yy, zz);
		//if (R < particles[i].calcR(Z/R)*0.999) return E;
	}
	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		T = acos(Z / R);
		P = atan2(Y, X);
		if (ml.nbrMediums == 1) {
			VSWF wf(N, particles[0].k1 * R, Z / R, P);
			ii = 0;
			for (int n = 1; n <= N; n++) {
				for (int m = -n; m <= n; m++) {
					wf.calcVSWFs(m, n, 3, true);
					E = E
							+ toCartesianCoordinates(
									Pmni.t[ii + 2 * NM * i] * wf.Mmn
											+ Pmni.t[ii + 2 * NM * i + NM]
													* wf.Nmn, T, P);
					ii++;
				}
			}
		} else {
			refl1.calcVSWFsIntegral(X, Y, zz, z[i]);
			int layJ = ml.getMedium(z[i]);
			int layI = ml.getMedium(zz);

			if (layI == layJ) {
				VSWF wf(N, ml.ki.t[layI] * R, Z / R, P);
				ii = 0;
				for (int n = 1; n <= N; n++) {
					for (int m = -n; m <= n; m++) {
						wf.calcVSWFs(m, n, 3, true);
						for (int iii = 0; iii < 3; iii++) {
							Mmn.t[iii] = refl1.Mmn.t[ii][iii];
							Nmn.t[iii] = refl1.Nmn.t[ii][iii];
						}
						E = E
								+ toCartesianCoordinates(
										Pmni.t[ii + 2 * NM * i] * wf.Mmn
												+ Pmni.t[ii + 2 * NM * i + NM]
														* wf.Nmn, T, P);
						E = E + Pmni.t[ii + 2 * NM * i] * Mmn
								+ Pmni.t[ii + 2 * NM * i + NM] * Nmn;
						ii++;
					}
				}
			} else {
				ii = 0;
				for (int n = 1; n <= N; n++) {
					for (int m = -n; m <= n; m++) {
						for (int iii = 0; iii < 3; iii++) {
							Mmn.t[iii] = refl1.Mmn.t[ii][iii];
							Nmn.t[iii] = refl1.Nmn.t[ii][iii];
						}
						E = E
								+ (Pmni.t[ii + 2 * NM * i] * Mmn
										+ Pmni.t[ii + 2 * NM * i + NM] * Nmn);
						ii++;
					}
				}
			}
		}
	}
	return E;
}

Vector3d MultipleScattering::calcHsca(const double &xx, const double &yy,
		const double &zz) {
	Vector3d E(0.);
	double X, Y, Z, R, T, P;
	int ii;
	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		//if (R < particles[i].calcR(Z/R)*0.9999999) return calcHint(i,xx,yy,zz);
		if (R < particles[i].calcR(Z / R) * 0.9999999)
			return E;
	}
	//return E;
	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		T = acos(Z / R);
		P = atan2(Y, X);
		VSWF wf(N, particles[0].k1 * R, Z / R, P);
		ii = 0;
		for (int n = 1; n <= N; n++) {
			for (int m = -n; m <= n; m++) {
				wf.calcVSWFs(m, n, 3, true);
				E = E
						+ toCartesianCoordinates(
								Pmni.t[ii + 2 * NM * i] * wf.Nmn
										+ Pmni.t[ii + 2 * NM * i + NM] * wf.Mmn,
								T, P);
				ii++;
			}
		}
	}
	return -i_ * sqrt(eps0 * particles[0].n1 * particles[0].n1 / mu0) * E;
}

Vector3d MultipleScattering::calcEscaInf(const double &theta,
		const double &phi) {
	Vector3d E(0.);
	int ii = 0;
	Complex cst;
	Complex dephasing;

	if (ml.nbrMediums == 1) {
		VSWFinf wf(N, theta, phi);
		for (int i = 0; i < Np; i++) {
			dephasing = exp(
					-i_ * particles[0].k1
							* (x[i] * cos(phi) * sin(theta)
									+ y[i] * sin(phi) * sin(theta)
									+ z[i] * cos(theta)));
			ii = 0;
			for (int n = 1; n <= N; n++) {
				cst = powi(-n - 1);
				for (int m = -n; m <= n; m++) {
					wf.calcVSWFs(m, n);
					E = E
							+ dephasing * cst
									* (Pmni.t[ii + 2 * i * NM] * wf.Mmn
											+ i_ * Pmni.t[ii + NM + 2 * i * NM]
													* wf.Nmn);
					ii++;
				}
			}
		}
		E.t[0] /= particles[0].k1;
		E.t[1] /= particles[0].k1;
		E.t[2] /= particles[0].k1;
	} else {
		Complex cstSnell;
		Vector app, asp, apt, ast, a(2, 1.);
		double zz;
		int lay;
		for (int i = 0; i < Np; i++) {
			lay = ml.getMedium(z[i]);
			if (theta <= piDblby2) {
				dephasing = exp(
						-i_ * ml.ki.t[ml.nbrMediums - 1]
								* (x[i] * cos(phi) * sin(theta)
										+ y[i] * sin(phi) * sin(theta)));
				zz = ml.Zi.t[ml.nbrMediums - 2] + 0.01e-9;
				ml.setKxy(
						ml.ki.t[ml.nbrMediums - 1]
								* sqrt(
										cos(phi) * sin(theta) * cos(phi)
												* sin(theta)
												+ sin(phi) * sin(theta)
														* sin(phi)
														* sin(theta)));

				VSWFinf wf(N, theta, phi);
				cstSnell = ml.ki.t[ml.nbrMediums - 1]
						* ml.calcKz(ml.nbrMediums - 1, true)
						/ (ml.ki.t[lay] * ml.calcKz(lay, true));
				ii = 0;
				for (int n = 1; n <= N; n++) {
					cst = dephasing * powi(-n - 1);
					for (int m = -n; m <= n; m++) {
						a.t[1] = powm1(n - m);
						app = ml.sourceInt(1, a, z[i], zz);
						asp = ml.sourceInt(0, a, z[i], zz);
						a.t[1] = -powm1(n - m);
						apt = ml.sourceInt(1, a, z[i], zz);
						ast = ml.sourceInt(0, a, z[i], zz);

						wf.calcVSWFs(m, n);
						E.t[1] += cst * cstSnell
								* (Pmni.t[ii + 2 * i * NM] * app.t[0]
										* wf.Mmn.t[1]
										+ i_ * Pmni.t[ii + NM + 2 * i * NM]
												* apt.t[0] * wf.Nmn.t[1]);
						E.t[2] += cst * cstSnell
								* (Pmni.t[ii + 2 * i * NM] * ast.t[0]
										* wf.Mmn.t[2]
										+ i_ * Pmni.t[ii + NM + 2 * i * NM]
												* asp.t[0] * wf.Nmn.t[2]);
						if (z[i] > zz)
							E =
									E
											+ cst
													* (Pmni.t[ii + 2 * i * NM]
															* wf.Mmn
															+ i_
																	* Pmni.t[ii
																			+ NM
																			+ 2
																					* i
																					* NM]
																	* wf.Nmn);
						ii++;
					}
				}
				E.t[0] /= ml.ki.t[ml.nbrMediums - 1];
				E.t[1] /= ml.ki.t[ml.nbrMediums - 1];
				E.t[2] /= ml.ki.t[ml.nbrMediums - 1];
			} else {
				dephasing = exp(
						-i_ * ml.ki.t[0]
								* (x[i] * cos(phi) * sin(theta)
										+ y[i] * sin(phi) * sin(theta)));
				zz = ml.Zi.t[0] - 0.01e-9;
				ml.setKxy(
						ml.ki.t[0]
								* sqrt(
										cos(phi) * sin(theta) * cos(phi)
												* sin(theta)
												+ sin(phi) * sin(theta)
														* sin(phi)
														* sin(theta)));
				VSWFinf wf(N, theta, phi);
				VSWFinf wf1(N, piDbl - theta, phi);

				cstSnell = ml.ki.t[0] * ml.calcKz(0, true)
						/ (ml.ki.t[lay] * ml.calcKz(lay, true));
				ii = 0;
				for (int n = 1; n <= N; n++) {
					cst = dephasing * powi(-n - 1);
					for (int m = -n; m <= n; m++) {
						a.t[1] = powm1(n - m);
						app = ml.sourceInt(1, a, z[i], zz);
						asp = ml.sourceInt(0, a, z[i], zz);
						a.t[1] = -powm1(n - m);
						apt = ml.sourceInt(1, a, z[i], zz);
						ast = ml.sourceInt(0, a, z[i], zz);
						wf.calcVSWFs(m, n);
						wf1.calcVSWFs(m, n);
						E.t[1] += cst * cstSnell
								* (Pmni.t[ii + 2 * i * NM] * app.t[1]
										* wf1.Mmn.t[1]
										+ i_ * Pmni.t[ii + NM + 2 * i * NM]
												* apt.t[1] * wf1.Nmn.t[1]);
						E.t[2] += cst * cstSnell
								* (Pmni.t[ii + 2 * i * NM] * ast.t[1]
										* wf1.Mmn.t[2]
										+ i_ * Pmni.t[ii + NM + 2 * i * NM]
												* asp.t[1] * wf1.Nmn.t[2]);
						if (z[i] < zz)
							E =
									E
											+ cst
													* (Pmni.t[ii + 2 * i * NM]
															* wf.Mmn
															+ i_
																	* Pmni.t[ii
																			+ NM
																			+ 2
																					* i
																					* NM]
																	* wf.Nmn);
						ii++;
					}
				}
				E.t[0] /= ml.ki.t[0];
				E.t[1] /= ml.ki.t[0];
				E.t[2] /= ml.ki.t[0];
			}
		}
	}
	return E;
}

Vector3d MultipleScattering::calcHscaInf(const double &theta,
		const double &phi) {
	Vector3d E(0.);
	VSWFinf wf(N, theta, phi);
	int ii = 0;
	for (int n = 1; n <= N; n++) {
		for (int m = -n; m <= n; m++) {
			wf.calcVSWFs(m, n);
			E = E + Pmn.t[ii] * wf.Nmn + Pmn.t[ii + NM] * wf.Mmn;
			ii++;
		}
	}
	E.t[0] /= particles[0].k1;
	E.t[1] /= particles[0].k1;
	E.t[2] /= particles[0].k1;
	return E;
}

Vector3d MultipleScattering::calcEint(const int &part, const double &xx,
		const double &yy, const double &zz) {
	Vector3d E(0.);
	double X, Y, Z, R, T, P;
	int ii;
	X = xx - x[part];
	Y = yy - y[part];
	Z = zz - z[part];
	R = sqrt(X * X + Y * Y + Z * Z);
	T = acos(Z / R);
	P = atan2(Y, X);
	VSWF wf(N, particles[part].k2 * R, Z / R, P);
	ii = 0;
	for (int n = 1; n <= N; n++) {
		for (int m = -n; m <= n; m++) {
			wf.calcVSWFs(m, n, 1, true);
			E = E
					+ toCartesianCoordinates(
							Cmni.t[ii + 2 * NM * part] * wf.Mmn
									+ Cmni.t[ii + 2 * NM * part + NM] * wf.Nmn,
							T, P);
			ii++;
		}
	}
	return E;
}

Vector3d MultipleScattering::calcHint(const int &part, const double &xx,
		const double &yy, const double &zz) {
	Vector3d E(0.);
	double X, Y, Z, R, T, P;
	int ii;
	X = xx - x[part];
	Y = yy - y[part];
	Z = zz - z[part];
	R = sqrt(X * X + Y * Y + Z * Z);
	T = acos(Z / R);
	P = atan2(Y, X);
	VSWF wf(N, particles[part].k2 * R, Z / R, P);
	ii = 0;
	for (int n = 1; n <= N; n++) {
		for (int m = -n; m <= n; m++) {
			wf.calcVSWFs(m, n, 1, true);
			E = E
					+ toCartesianCoordinates(
							Cmni.t[ii + 2 * NM * part] * wf.Nmn
									+ Cmni.t[ii + 2 * NM * part + NM] * wf.Mmn,
							T, P);
			ii++;
		}
	}
	return -i_ * sqrt(eps0 * particles[part].n2 * particles[part].n2 / mu0) * E;
}

Vector3d MultipleScattering::calcEinc(const double &xx, const double &yy,
		const double &zz) {
	Vector3d E(0.), pos;
	double X, Y, Z, R;

	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		if (R < particles[i].calcR(Z / R) * 0.9999999)
			return E;
	}
	int lay = ml.getMedium(zz);
	pos.t[0] = xx;
	pos.t[1] = yy;
	pos.t[2] = zz;
	incidence.setPlaneWaveCoeffs(N, pos);

	Complex bP = acos(incidence.ml.calcKz(lay, true) / incidence.ml.ki.t[lay]),
			bM = piDbl
					- acos(
							incidence.ml.calcKz(lay, true)
									/ incidence.ml.ki.t[lay]);
	Vector3d EP(0.), EM(0.);
	Vector ap = incidence.ml.sourceExtIn(incidence.aincp, 1, zz);
	Vector as = incidence.ml.sourceExtIn(incidence.aincs, 0, zz);
	ap = exp(
			i_ * incidence.ml.kxy
					* (xx * cos(incidence.alpha0) + yy * sin(incidence.alpha0)))
			* ap;
	as = exp(
			i_ * incidence.ml.kxy
					* (xx * cos(incidence.alpha0) + yy * sin(incidence.alpha0)))
			* as;
	EP.t[1] = ap.t[0];
	EP.t[2] = as.t[0];
	EM.t[1] = ap.t[1];
	EM.t[2] = as.t[1];
	//cout << incidence.aincp << endl;
	EP = toCartesianCoordinates(EP, bP, incidence.alpha0);
	EM = toCartesianCoordinates(EM, bM, incidence.alpha0);

	return EM + EP;

	//E.t[1] = incidence.Eb*exp(i_*particles[0].k1*((xx*cos(incidence.alpha0)+yy*sin(incidence.alpha0))*sin(incidence.beta0)+zz*cos(incidence.beta0)));
	//E.t[2] = incidence.Ea*exp(i_*particles[0].k1*((xx*cos(incidence.alpha0)+yy*sin(incidence.alpha0))*sin(incidence.beta0)+zz*cos(incidence.beta0)));
	//E = toCartesianCoordinates(E,incidence.beta0,incidence.alpha0);
	//return E;
}

Vector3d MultipleScattering::calcHinc(const double &xx, const double &yy,
		const double &zz) {
	Vector3d E(0.);
	double X, Y, Z, R;
	for (int i = 0; i < Np; i++) {
		X = xx - x[i];
		Y = yy - y[i];
		Z = zz - z[i];
		R = sqrt(X * X + Y * Y + Z * Z);
		if (R < particles[i].calcR(Z / R) * 0.9999999)
			return E;
	}
	E.t[2] = incidence.Eb
			* exp(
					i_ * particles[0].k1
							* ((xx * cos(incidence.alpha0)
									+ yy * sin(incidence.alpha0))
									* sin(incidence.beta0)
									+ zz * cos(incidence.beta0)));
	E.t[1] = incidence.Ea
			* exp(
					i_ * particles[0].k1
							* ((xx * cos(incidence.alpha0)
									+ yy * sin(incidence.alpha0))
									* sin(incidence.beta0)
									+ zz * cos(incidence.beta0)));
	E = toCartesianCoordinates(E, incidence.beta0, incidence.alpha0);
	return -i_ * sqrt(eps0 * particles[0].n1 * particles[0].n1 / mu0) * E;
}

Vector MultipleScattering::calcStokesParametersInc() {
	Vector I(8);

	Complex cst1 = 1.;
	Complex cst2 = 1.;
	if (incidence.beta0 < piDbl / 2.)
		cst1 *= incidence.ml.calcKz(ml.nbrMediums - 1, true)
				/ incidence.ml.calcKz(0, true);
	else
		cst2 *= incidence.ml.calcKz(0, true)
				/ incidence.ml.calcKz(ml.nbrMediums - 1, true);
	I.t[0] = cst1 * (norm(incidence.aoutp.t[1]) + norm(incidence.aouts.t[1]));
	I.t[1] = cst1 * (norm(incidence.aoutp.t[1]) - norm(incidence.aouts.t[1]));
	I.t[2] = cst1
			* (-incidence.aouts.t[1] * conj(incidence.aoutp.t[1])
					- incidence.aoutp.t[1] * conj(incidence.aouts.t[1]));
	I.t[3] = cst1 * i_
			* (incidence.aouts.t[1] * conj(incidence.aoutp.t[1])
					- incidence.aoutp.t[1] * conj(incidence.aouts.t[1]));

	I.t[4] = cst2 * (norm(incidence.aoutp.t[0]) + norm(incidence.aouts.t[0]));
	I.t[5] = cst2 * (norm(incidence.aoutp.t[0]) - norm(incidence.aouts.t[0]));
	I.t[6] = cst2
			* (-incidence.aouts.t[0] * conj(incidence.aoutp.t[0])
					- incidence.aoutp.t[0] * conj(incidence.aouts.t[0]));
	I.t[7] = cst2 * i_
			* (incidence.aouts.t[0] * conj(incidence.aoutp.t[0])
					- incidence.aoutp.t[0] * conj(incidence.aouts.t[0]));
	return I;
}

Vector MultipleScattering::calcStokesParametersExt() {
	Vector I(8);
	double ctp = real(
			incidence.ml.calcKz(ml.nbrMediums - 1, true)
					/ ml.ki.t[ml.nbrMediums - 1]);
	double ctm = real(incidence.ml.calcKz(0, false) / ml.ki.t[0]);
	//double ctp,ctm;
	double alpha = incidence.alpha0;
	Complex cstp, cstm;
	cstp = piDbl / ml.ki.t[ml.nbrMediums - 1];
	cstm = piDbl / ml.ki.t[0];
	Vector3d Ep, Em;
	Ep = calcEscaInf(acos(ctp), alpha);
	Em = calcEscaInf(acos(ctm), alpha);
	if (incidence.beta0 < piDbl / 2.)
		cstp *= incidence.ml.calcKz(ml.nbrMediums - 1, true)
				/ incidence.ml.calcKz(0, true);
	else
		cstm *= incidence.ml.calcKz(0, true)
				/ incidence.ml.calcKz(ml.nbrMediums - 1, true);
	I.t[0] = cstp
			* imag(
					Ep.t[1] * conj(incidence.aoutp.t[1])
							+ Ep.t[2] * conj(incidence.aouts.t[1]));
	I.t[1] = cstp
			* imag(
					Ep.t[1] * conj(incidence.aoutp.t[1])
							- Ep.t[2] * conj(incidence.aouts.t[1]));
	I.t[2] = -cstp
			* imag(
					Ep.t[1] * conj(incidence.aouts.t[1])
							+ Ep.t[2] * conj(incidence.aoutp.t[1]));
	I.t[3] = -cstp
			* real(
					Ep.t[1] * conj(incidence.aouts.t[1])
							- Ep.t[2] * conj(incidence.aoutp.t[1]));

	I.t[4] = cstm
			* imag(
					Em.t[1] * conj(incidence.aoutp.t[0])
							+ Em.t[2] * conj(incidence.aouts.t[0]));
	I.t[5] = cstm
			* imag(
					Em.t[1] * conj(incidence.aoutp.t[0])
							- Em.t[2] * conj(incidence.aouts.t[0]));
	I.t[6] = -cstm
			* imag(
					Em.t[1] * conj(incidence.aouts.t[0])
							+ Em.t[2] * conj(incidence.aoutp.t[0]));
	I.t[7] = -cstm
			* real(
					Em.t[1] * conj(incidence.aouts.t[0])
							- Em.t[2] * conj(incidence.aoutp.t[0]));
	return I;
}

Matrix MultipleScattering::calcScatteringMatrix(const double &ctp,
		const double &ctm, const double &alphap, const double &alpham) {
	Matrix S(2, 4, 0.);
	Vector3d Ep, Em;
	Complex Ea = incidence.Ea, Eb = incidence.Eb;
	//if (Eb != 0.) {
	incidence.Ea = 0.;
	incidence.Eb = 1.;
	solve();

	Ep = calcEscaInf(acos(ctp), alphap);
	Em = calcEscaInf(acos(ctm), alpham);
	//cout << Ep << "\n" << Em << endl;
	S.t[0][0] = Ep.t[1];
	S.t[1][0] = Ep.t[2];
	S.t[0][2] = Em.t[1];
	S.t[1][2] = Em.t[2];
	//}

	incidence.Ea = 1.;
	incidence.Eb = 0.;
	solve();

	Ep = calcEscaInf(acos(ctp), alphap);
	Em = calcEscaInf(acos(ctm), alpham);
	S.t[0][1] = Ep.t[1];
	S.t[1][1] = Ep.t[2];
	S.t[0][3] = Em.t[1];
	S.t[1][3] = Em.t[2];
	incidence.Ea = Ea;
	incidence.Eb = Eb;

	//cout << S << "\n\n" << endl;
	return S;
}

Matrix MultipleScattering::calcExtinctionMatrix() {
	Vector3d Ep, Em;
	Matrix S = calcScatteringMatrix(
			real(
					incidence.ml.calcKz(ml.nbrMediums - 1, true)
							/ ml.ki.t[ml.nbrMediums - 1]),
			real(incidence.ml.calcKz(0, false) / ml.ki.t[0]), incidence.alpha0,
			incidence.alpha0);
	Matrix K(4, 8);
	Complex cstp = 2. * piDbl / ml.ki.t[ml.nbrMediums - 1], cstm = 2. * piDbl
			/ ml.ki.t[0];
	for (int i = 0; i < 4; i++) {
		K.t[i][i] = cstp * imag(S.t[0][0] + S.t[1][1]);
		K.t[i][i + 4] = cstm * imag(S.t[0][2] + S.t[1][3]);
	}
	K.t[0][1] = K.t[1][0] = cstp * imag(S.t[0][0] - S.t[1][1]);
	K.t[0][5] = K.t[1][4] = cstm * imag(S.t[0][2] - S.t[1][3]);

	K.t[0][2] = K.t[2][0] = -cstp * imag(S.t[0][1] + S.t[1][0]);
	K.t[0][6] = K.t[2][4] = -cstm * imag(S.t[0][3] - S.t[1][2]);

	K.t[0][3] = K.t[3][0] = cstp * real(S.t[1][0] - S.t[0][1]);
	K.t[0][7] = K.t[3][4] = cstm * real(S.t[1][2] - S.t[0][3]);

	K.t[1][2] = cstp * imag(S.t[1][0] - S.t[0][1]);
	K.t[2][1] = -cstp * imag(S.t[1][0] - S.t[0][1]);
	K.t[1][6] = cstm * imag(S.t[1][2] - S.t[0][3]);
	K.t[2][5] = -cstm * imag(S.t[1][2] - S.t[0][3]);

	K.t[1][3] = -cstp * real(S.t[0][1] + S.t[1][0]);
	K.t[3][1] = cstp * real(S.t[0][1] + S.t[1][0]);
	K.t[1][7] = -cstm * real(S.t[0][3] + S.t[1][2]);
	K.t[3][5] = cstm * real(S.t[0][3] + S.t[1][2]);

	K.t[2][3] = cstp * real(S.t[1][1] - S.t[0][0]);
	K.t[3][2] = -cstp * real(S.t[1][1] - S.t[0][0]);
	K.t[2][7] = cstm * real(S.t[1][3] - S.t[0][2]);
	K.t[3][6] = -cstm * real(S.t[1][3] - S.t[0][2]);

	return K;
}
