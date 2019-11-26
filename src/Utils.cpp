#include "Utils.h"

int ii(const int &n, const int &m) {
	return n * n + n + m - 1;
}

double fact(const int &n) {
	if (n == 0 || n == 1)
		return 1.;
	double x = 1.;
	for (int i = 2; i <= n; i++)
		x *= double(i);
	return x;
}

double dblfact(const int &n) {
	if (n <= 1)
		return 1.;
	double x = double(n);
	for (int i = n - 2; i >= 1; i -= 2)
		x *= double(i);
	return x;
}

double binomialCoef(const int &n, const int &l) {
	int nl = max(n, l);
	double lns[nl + 1];
	lns[0] = lns[1] = 0.;
	for (int i = 1; i < nl; i++)
		lns[i + 1] = lns[i] + log(double(i + 1));
	return exp(lns[n] - lns[l] - lns[n - l]);
}

double powm1(const int &n) {
	if ((n % 2) == 0)
		return 1.;
	return -1.;
}

void swap(int &a, int &b) {
	int c = a;
	a = b;
	b = c;
}

void swap(double &a, double &b) {
	double c = a;
	a = b;
	b = c;
}

double max(const double &a, const double &b) {
	if (a > b)
		return a;
	return b;
}

int max(const int &a, const int &b) {
	if (a > b)
		return a;
	return b;
}
