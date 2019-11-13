#ifndef SPECIALFUNCTIONS_H_INCLUDED
#define SPECIALFUNCTIONS_H_INCLUDED

#include "Complex.h"

double gammaln(const double &xx);
Complex fresnelCoefRP(const Complex &m,const Complex &cosBeta1);
Complex fresnelCoefRS(const Complex &m,const Complex &cosbeta1);

Complex fresnelCoefTP(const Complex &m,const Complex &cosBeta1);
Complex fresnelCoefTS(const Complex &m,const Complex &cosbeta1);

double gegenbauer(const double &alph, const int &n, const double &x);

Complex Amnplus(const int &m, const int &n, const double &v);
Complex Bmnplus(const int &m, const int &n, const double &v);
#endif // SPECIALFUNCTIONS_H_INCLUDED
