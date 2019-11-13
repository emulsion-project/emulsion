#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include "RealContainers.h"
#include "Matrix.h"
#include "SpecialFunctions.h"

RMatrix gaussLegQuad(const int &n,const double &x1, const double &x2);
Matrix gaussLegQuad(const int &n,const Complex &x1, const Complex &x2);
RMatrix gaussLagQuad(const int &n);
RMatrix gaussLagQuad(const int &n,const double &alf);

#endif // INTEGRATION_H_INCLUDED
