#ifndef LEGENDREFUNCTIONS_H_INCLUDED
#define LEGENDREFUNCTIONS_H_INCLUDED

#include "Vector.h"
#include "RealContainers.h"

class LegendreFunctionsC {
public:
	LegendreFunctionsC();
	LegendreFunctionsC(const int &NN, const Complex &xx);
	void assign(const int &NN, const Complex &xx);
	void calcElmts();
	void calcElmts0();
	void calcElmtsPi();
	~LegendreFunctionsC();
	int N;
	Complex x;
	Vector *pmn;
	Vector *pimn;
	Vector *taumn;
	Complex Pmn(const int &n, const int &m);
	Complex Pimn(const int &n, const int &m);
	Complex Taumn(const int &n, const int &m);
};

class LegendreFunctions {
public:
	LegendreFunctions();
	LegendreFunctions(const int &NN, const double &xx);
	void assign(const int &NN, const double &xx);
	void calcElmts();
	void calcElmts0();
	void calcElmtsPi();
	~LegendreFunctions();
	int N;
	double x;
	RVector *pmn;
	RVector *pimn;
	RVector *taumn;
	double Pmn(const int &n, const int &m);
	double Pimn(const int &n, const int &m);
	double Taumn(const int &n, const int &m);
};

#endif // LEGENDREFUNCTIONS_H_INCLUDED
