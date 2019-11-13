#ifndef TRANSMISSION_H_INCLUDED
#define TRANSMISSION_H_INCLUDED

#include "Matrix.h"
#include "BesselFunctions.h"
#include "LegendreFunctions.h"
#include "Integration.h"
#include "Multilayer.h"

class Transmission
{
public:
  Transmission();
  Transmission(const double &lambd, const int &n, const int &ng, const double &xx, const double &yy, const double &zzi, const double zzj, const bool &self);
  void setParameters(const double &lambd, const int &n, const int &ng, const double &xx, const double &yy, const double &zzi, const double zzj, const bool &self);
  Matrix calcR(const bool &Up);
  Matrix calcRSelf(const bool &Up);
  bool self;
  int N,NM,Ng,lay;
  double dx,dy,zi,zj,xMax,pasX;
  Complex *xInt;
  Matrix R;
  Multilayer ml;
  bool Self;
  RMatrix gaussLag;
  RMatrix gaussLeg;
  LegendreFunctions *leg;
  LegendreFunctionsC *legC;
  Bessel bes;
};

#endif // TRANSMISSION_H_INCLUDED
