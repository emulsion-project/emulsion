#include "ClebshGordan.h"

ClebshGordan::ClebshGordan() {
	N = 0;
	coefs = NULL;
}

ClebshGordan::~ClebshGordan() {
	if (coefs != NULL) {
		for (int i = 0; i < (N + 1) * (N + 1); i++)
			delete[] (coefs[i]);
		delete[] (coefs);
	}
}

ClebshGordan::ClebshGordan(const int &NN) {
	N = NN;
	const int R = (N + 1) * (N + 1);
	int row, col, n, l, m, k, w, lim;
	double a, b, c;
	coefs = new RVector*[R];
	for (int r = 0; r < R; r++)
		coefs[r] = new RVector[R];
	row = col = 0;
	for (n = 0; n <= N; n++) {
		for (m = -n; m <= n; m++) {
			col = 0;
			for (l = 0; l <= N; l++) {
				for (k = -l; k <= l; k++) {
					lim = max(abs(n - l), abs(m + k));
					coefs[row][col] = RVector(n + l - lim + 1, 0.);
					coefs[row][col].t[n + l - lim] = sqrt(
							binomialCoef(n + l + m + k, l + k)
									* binomialCoef(n + l - m - k, l - k)
									/ binomialCoef(n + n + l + l, l + l));
					w = n + l + 1;
					if (w - 2 >= lim) {
						a = double(w + w - 2)
								* sqrt(
										double(w + w - 3) * double(w + w - 1)
												/ (double(w - m - k - 1)
														* double(w + m + k - 1)
														* double(l - n + w - 1)
														* double(n - l + w - 1)
														* double(n + l - w + 2)
														* double(n + l + w)));
						b = (double(m - k) * double(w * w - w)
								- double(m + k)
										* (double(n * n + n - l * l - l)))
								/ double(2 * w * w - 2 * w);
						coefs[row][col].t[w - 2 - lim] = a * b
								* coefs[row][col].t[w - 1 - lim];
					}

					for (w = n + l; w >= lim + 2; w--) {
						a = double(w + w - 2)
								* sqrt(
										double(w + w - 3) * double(w + w - 1)
												/ (double(w - m - k - 1)
														* double(w + m + k - 1)
														* double(l - n + w - 1)
														* double(n - l + w - 1)
														* double(n + l - w + 2)
														* double(n + l + w)));
						b = (double(m - k) * double(w * w - w)
								- double(m + k)
										* (double(n * n + n - l * l - l)))
								/ double(2 * w * w - 2 * w);
						c =
								sqrt(
										double(w + m + k) * double(w - m - k)
												* double(n - l + w)
												* double(l - n + w)
												* double(n + l - w + 1)
												* double(n + l + w + 1)
												/ (double(w + w + 1)
														* double(w + w - 1)))
										/ double(w + w);
						coefs[row][col].t[w - 2 - lim] = a
								* (b * coefs[row][col].t[w - 1 - lim]
										- c * coefs[row][col].t[w - lim]);
					}
					col++;
				}
			}
			row++;
		}
	}
}

void ClebshGordan::assign(const int &NN) {
	if (N == NN)
		return;
	if (coefs != NULL) {
		for (int i = 0; i < (N + 1) * (N + 1); i++)
			delete[] (coefs[i]);
		delete[] (coefs);
	}
	N = NN;
	const int R = (N + 1) * (N + 1);
	int row, col, n, l, m, k, w, lim;
	double a, b, c;
	coefs = new RVector*[R];
	for (int r = 0; r < R; r++)
		coefs[r] = new RVector[R];
	row = col = 0;
	for (n = 0; n <= N; n++) {
		for (m = -n; m <= n; m++) {
			col = 0;
			for (l = 0; l <= N; l++) {
				for (k = -l; k <= l; k++) {
					lim = max(abs(n - l), abs(m + k));
					coefs[row][col] = RVector(n + l - lim + 1, 0.);
					coefs[row][col].t[n + l - lim] = sqrt(
							binomialCoef(n + l + m + k, l + k)
									* binomialCoef(n + l - m - k, l - k)
									/ binomialCoef(n + n + l + l, l + l));
					w = n + l + 1;
					if (w - 2 >= lim) {
						a = double(w + w - 2)
								* sqrt(
										double(w + w - 3) * double(w + w - 1)
												/ (double(w - m - k - 1)
														* double(w + m + k - 1)
														* double(l - n + w - 1)
														* double(n - l + w - 1)
														* double(n + l - w + 2)
														* double(n + l + w)));
						b = (double(m - k) * double(w * w - w)
								- double(m + k)
										* (double(n * n + n - l * l - l)))
								/ double(2 * w * w - 2 * w);
						coefs[row][col].t[n + l - 1 - lim] = a * b
								* coefs[row][col].t[n + l - lim];
					}
					for (w = n + l; w >= lim + 2; w--) {
						a = double(w + w - 2)
								* sqrt(
										double(w + w - 3) * double(w + w - 1)
												/ (double(w - m - k - 1)
														* double(w + m + k - 1)
														* double(l - n + w - 1)
														* double(n - l + w - 1)
														* double(n + l - w + 2)
														* double(n + l + w)));
						b = (double(m - k) * double(w * w - w)
								- double(m + k)
										* (double(n * n + n - l * l - l)))
								/ double(2 * w * w - 2 * w);
						c = 1. / double(w + w)
								* sqrt(
										double(w + m + k) * double(w - m - k)
												* double(n - l + w)
												* double(l - n + w)
												* double(n + l - w + 1)
												* double(n + l + w + 1)
												/ (double(w + w + 1)
														* double(w + w - 1)));
						coefs[row][col].t[w - 2 - lim] = a
								* (b * coefs[row][col].t[w - 1 - lim]
										- c * coefs[row][col].t[w - lim]);
					}
					col++;
				}
			}
			row++;
		}
	}
}

int ClebshGordan::mn(const int &n, const int &m) {
	return n * n + n + m;
}

double ClebshGordan::coef(const int &w, const int &n, const int &m,
		const int &l, const int &k) {
	if ((abs(m) > n) || (abs(k) > l) || (abs(m + k) > w))
		return 0.;
	if ((abs(n - l) > w) || (w > n + l))
		return 0.;
	int lim = max(abs(n - l), abs(m + k));
	return coefs[mn(n, m)][mn(l, k)].t[w - lim];
}
