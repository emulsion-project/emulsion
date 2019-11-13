#ifndef REFRACTIVEINDEX_H_INCLUDED
#define REFRACTIVEINDEX_H_INCLUDED

#include "Complex.h"
#include <fstream>

class RefractiveIndex
{
public:
  int type,material;
  Complex n;
  RefractiveIndex();
  RefractiveIndex(const int &mat, const double &R);
  void setMaterial(const int &mat, const double &R);
  string fileName,khiInterFile;
  double wP,gamma0,A,vF,R,lambda;
  Complex calcIndex(const double &lambd);
  Complex interpolation(void);
  Complex drude();
  Complex drudeMod();
};

#endif // REFRACTIVEINDEX_H_INCLUDED
