#ifndef WIGNERFUNCTIONS_H_INCLUDED
#define WIGNERFUNCTIONS_H_INCLUDED

#include "LegendreFunctions.h"

class WignerFunctions
{
public:
  WignerFunctions();
  WignerFunctions(const int &n, const double &thet);
  ~WignerFunctions();
  void assign(const int &nn, const double &thet);
  double dlmn(const int &l, const int &m, const int &n);
  double x,theta;
  int N,n0;
  RVector **coefs;
};

class WignerFunctionsC
{
public:
  WignerFunctionsC();
  WignerFunctionsC(const int &n, const Complex &thet);
  ~WignerFunctionsC();
  void assign(const int &nn, const Complex &thet);
  Complex dlmn(const int &l, const int &m, const int &n);
  Complex x,theta;
  int N,n0;
  Vector **coefs;
};

#endif // WIGNERFUNCTIONS_H_INCLUDED
