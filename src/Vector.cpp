#include "Vector.h"
#include <iostream>
using namespace std;

/*********************************
 *     Class constructors        *
 *********************************/

Vector::Vector() :
		n(0), t(NULL) {
}
Vector::Vector(const int &n) :
		n(n), t(n > 0 ? new Complex[n] : NULL) {
}
Vector::Vector(const int &n, const Complex &x) :
		n(n), t(n > 0 ? new Complex[n] : NULL) {
	for (int i = 0; i < n; i++)
		t[i] = x;
}
Vector::Vector(const Vector &v) :
		n(v.n), t(n > 0 ? new Complex[n] : NULL) {
	for (int i = 0; i < n; i++)
		t[i] = v.t[i];
}

/*********************************
 *     Class destructor          *
 *********************************/

Vector::~Vector() {
	if (t != NULL)
		delete[] (t);
}

/*********************************
 *     Class member operators    *
 *********************************/

Vector& Vector::operator=(const Vector &v) {
	if (this != &v) {
		if (n != v.n) {
			if (t != NULL)
				delete[] (t);
			n = v.n;
			t = n > 0 ? new Complex[n] : NULL;
		}
		for (int i = 0; i < n; i++)
			t[i] = v.t[i];
	}
	return *this;
}

/*********************************
 *     Class member functions    *
 *********************************/

void Vector::resize(const int &m) {
	if (m != n) {
		if (t != NULL)
			delete[] (t);
		n = m;
		t = n > 0 ? new Complex[n] : NULL;
	}
}

void Vector::assign(const int &m, const Complex &z) {
	if (m != n) {
		if (t != NULL)
			delete[] (t);
		n = m;
		t = n > 0 ? new Complex[n] : NULL;
	}
	for (int i = 0; i < n; i++)
		t[i] = z;
}

/*********************************
 *     Associated operators      *
 *********************************/

Vector operator+(const Vector &v1, const Vector &v2) {
	Vector v3(v1.n);
	if (v1.n == v2.n)
		for (int i = 0; i < v3.n; i++)
			v3.t[i] = v1.t[i] + v2.t[i];
	return v3;
}

Vector operator-(const Vector &v1, const Vector &v2) {
	Vector v3(v1.n);
	if (v1.n == v2.n)
		for (int i = 0; i < v3.n; i++)
			v3.t[i] = v1.t[i] - v2.t[i];
	return v3;
}

Complex operator*(const Vector &v1, const Vector &v2) {
	Complex t = 0.;
	if (v1.n == v2.n)
		for (int i = 0; i < v1.n; i++)
			t += conj(v1.t[i]) * v2.t[i];
	return t;
}

Vector operator*(const Vector &v1, const Complex &z) {
	Vector v2(v1.n);
	for (int i = 0; i < v2.n; i++)
		v2.t[i] = v1.t[i] * z;
	return v2;
}

Vector operator*(const Complex &z, const Vector &v1) {
	Vector v2(v1.n);
	for (int i = 0; i < v2.n; i++)
		v2.t[i] = v1.t[i] * z;
	return v2;
}

std::ostream& operator<<(std::ostream &out, const Vector &v) {
	for (int i = 0; i < v.n; i++)
		out << v.t[i] << "\n";
	return out;
}

/*********************************
 *     Associated functions      *
 *********************************/

Complex max(const Vector &V) {
	Complex v = V.t[0];
	for (int i = 0; i < V.n; i++)
		if (v < V.t[i])
			v = V.t[i];
	return v;
}

Complex min(const Vector &V) {
	Complex v = V.t[0];
	for (int i = 0; i < V.n; i++)
		if (v > V.t[i])
			v = V.t[i];
	return v;
}

Complex sum(const Vector &V) {
	Complex z = 0.;
	for (int i = 0; i < V.n; i++)
		z += V.t[i];
	return z;
}

Vector poleAmplitude(const Vector &x, const Vector &y) {
	int n = x.n;
	Vector xpap(2, 0.);
	Complex aa = 0., bb = 0., a, b;
	for (int i = 0; i < n; i++) {
		a = x.t[i] * y.t[i];
		b = y.t[i];
		for (int j = 0; j < n; j++) {
			if (j != i) {
				a /= x.t[i] - x.t[j];
				b /= x.t[i] - x.t[j];
			}
		}
		aa += a;
		bb += b;
	}
	xpap.t[0] = aa / bb;
	aa = 0.;
	for (int i = 0; i < n; i++) {
		a = (x.t[i] - xpap.t[0]) * y.t[i];
		for (int j = 0; j < n; j++) {
			if (j != i)
				a *= (xpap.t[0] - x.t[j]) / (x.t[i] - x.t[j]);
		}
		aa += a;
	}
	xpap.t[1] = aa;
	return xpap;
}

double norm(const Vector &v) {
	double x = 0.;
	for (int i = 0; i < v.n; i++)
		x += norm(v.t[i]);
	return sqrt(x);
}

Vector fft(const Vector &v, const int &isign) {
	double *data;
	int nn, mmax, m, j, istep, i, n = v.n;
	double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
	Vector transf(n);

	data = new double[2 * n];
	for (i = 0; i < n; i++) {
		data[2 * i] = v.t[i].re;
		data[2 * i + 1] = v.t[i].im;
	}
	nn = n << 1;
	j = 1;
	for (i = 1; i < nn; i += 2) {
		if (j > i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (nn > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / double(mmax));
		wtemp = sin(0.5 * theta);
		wpr = -2. * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.;
		wi = 0.;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= nn; i += istep) {
				j = i + mmax;
				tempr = wr * data[j - 1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j - 1];
				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
	double mult = 1.;
	if (isign == -1)
		mult = 1. / double(n);
	transf.t[0] = Complex(data[0] * mult, data[1] * mult);
	for (i = 1; i < n; i++)
		transf.t[i] = mult
				* Complex(data[2 * n - 1 - 2 * (i - 1) - 1],
						data[2 * n - 1 - 2 * (i - 1)]);
	delete[] data;
	return transf;
}
