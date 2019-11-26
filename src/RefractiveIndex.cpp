#include "RefractiveIndex.h"
#include <iostream>
using namespace std;

RefractiveIndex::RefractiveIndex() {
}

RefractiveIndex::RefractiveIndex(const int &mat, const double &r) {
	ifstream ifs("Materials/materials.txt");
	material = mat;
	R = r;
	int i = 1;
	string str;
	while (i < material) {
		getline(ifs, str);
		i++;
	}
	ifs >> type;
	if (type == 0) {
		double N, K;
		ifs >> N >> K;
		n = Complex(N, K);
	} else if (type == 1)
		ifs >> fileName;
	else if (type == 2)
		ifs >> wP >> gamma0;
	else if (type == 3)
		ifs >> wP >> gamma0 >> A >> vF >> fileName;
}

void RefractiveIndex::setMaterial(const int &mat, const double &r) {
	ifstream ifs("Materials/materials.txt");
	material = mat;
	R = r;
	int i = 1;
	string str;
	while (i < material) {
		getline(ifs, str);
		i++;
	}
	ifs >> type;
	if (type == 0) {
		double N, K;
		ifs >> N >> K;
		n = Complex(N, K);
	} else if (type == 1)
		ifs >> fileName;
	else if (type == 2)
		ifs >> wP >> gamma0;
	else if (type == 3)
		ifs >> wP >> gamma0 >> A >> vF >> fileName;
	return;
}

Complex RefractiveIndex::calcIndex(const double &lambd) {
	lambda = lambd;
	if (type == 0)
		return n;
	else if (type == 1)
		return interpolation();
	else if (type == 2)
		return drude();
	return drudeMod();
}

Complex RefractiveIndex::interpolation() {
	ifstream ifs(fileName.c_str());
	int nbr;
	ifs >> nbr;
	int i = 0;
	double l1, l2, n1, n2, k1, k2, N, K;
	ifs >> l1 >> n1 >> k1;
	ifs >> l2 >> n2 >> k2;
	while (i < nbr - 1) {
		if ((l1 >= lambda && l2 <= lambda) || (l1 <= lambda && l2 >= lambda)) {
			N = (l2 - lambda) / (l2 - l1) * n1 + (lambda - l1) / (l2 - l1) * n2;
			K = (l2 - lambda) / (l2 - l1) * k1 + (lambda - l1) / (l2 - l1) * k2;
			break;
		}
		l1 = l2;
		n1 = n2;
		k1 = k2;
		ifs >> l2 >> n2 >> k2;
		i++;
	}
	return Complex(N, K);
}

Complex RefractiveIndex::drude() {
	double w = 2. * piDbl * c / lambda;
	double wwP = wP * e / hbar;
	return sqrt(-wwP * wwP / (w * w + i_ * w * gamma0 * e / hbar));
}

Complex RefractiveIndex::drudeMod() {
	double w = 2. * piDbl * c / lambda;
	double wwP = wP * e / hbar;
	double gamma = gamma0 * e / hbar + A * vF / R;
	Complex khiInter = interpolation();
	return sqrt(khiInter + 1. - wwP * wwP / (w * w + i_ * w * gamma));
}
