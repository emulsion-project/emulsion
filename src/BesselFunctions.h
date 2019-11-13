#ifndef BESSELFUNCTIONS_H_INCLUDED
#define BESSELFUNCTIONS_H_INCLUDED

#include "Vector.h"
#include "Integration.h"

class BesselFunctions
{
public:
  BesselFunctions();
  void calcCJ(const int &NN, const double &alph, const Complex &zz);
  void calcCJ(const int &NN, const Complex &zz);
  void calcJ(const int &NN, const double &alpha, const double &xx);
  void calcJ(const int &NN, const double &xx);
  void calcY(const bool &deriv);
  void calcH1(const bool &deriv);
  void calcH2(const bool &deriv);
  void calcSpherJ(const bool &deriv);
  void calcSpherY(const bool &deriv);
  void calcSpherH1(const bool &deriv);
  void calcSpherH2(const bool &deriv);
  Vector CJ,CdJ,CY,CdY,CH1,CdH1,CH2,CdH2;
  RVector J,dJ,Y,dY,H1,dH1,H2,dH2;
  int n,eta;
  Complex z;
  double x;
};


class SphericalBessel
{
public:
  SphericalBessel();
  SphericalBessel(const int &NN,const Complex &zz);
  void assign(const int &NN,const Complex &zz);
  Vector J,dJ,Y,dY;
  Complex jn(const int &n);
  Complex yn(const int &n);
  Complex h1n(const int &n);
  Complex h2n(const int &n);
  Complex djn(const int &n);
  Complex dyn(const int &n);
  Complex dh1n(const int &n);
  Complex dh2n(const int &n);
  int NP;
};

class Bessel
{
public:
  Bessel();
  Bessel(const int &NN,const Complex &zz);
  void assign(const int &NN, const Complex &zz);
  Vector J,Y,H1;
  RMatrix wx;
  Matrix wx1;
  int Ng;
  Complex Jn(const int &n);
};

Complex besselJn(const int &n, const Complex &z);
double besselK0(const double &x);
double besselK1(const double &x);
double besselKn(const int &n, const double &x);
double besselI0(const double &x);
double besselI1(const double &x);
double besselIn(const int &n, const double &x);

#endif // BESSELFUNCTIONS_H_INCLUDED
