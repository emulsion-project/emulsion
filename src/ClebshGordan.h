#ifndef CLEBSHGORDAN_H_INCLUDED
#define CLEBSHGORDAN_H_INCLUDED

#include "RealContainers.h"
#include "Utils.h"

class ClebshGordan {
public:
	ClebshGordan();
	ClebshGordan(const int &NN);
	~ClebshGordan();
	void assign(const int &NN);
	double calcCoef(const int &w, const int &n, const int &m, const int &l,
			const int &k);
	double coef(const int &w, const int &n, const int &m, const int &l,
			const int &k);
	int mn(const int &n, const int &m);
	int N;
	RVector **coefs;
};

#endif // CLEBSHGORDAN_H_INCLUDED
