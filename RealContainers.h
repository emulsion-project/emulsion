#ifndef REALCONTAINERS_H_INCLUDED
#define REALCONTAINERS_H_INCLUDED

#include <ostream>

/*********************************
 *     Class declaration         *
 *********************************/

class RVector
{
public:
  RVector();
  RVector(const int &n);
  RVector(const int &n,const double &x);
  RVector(const RVector &v);
  ~RVector();
  void resize(const int &m);
  void assign(const int &m,const double &x);
  RVector& operator=(const RVector &v);
  int n;
  double *t;
};

class IVector
{
public:
  IVector();
  IVector(const int &n);
  IVector(const int &n,const int &x);
  IVector(const IVector &v);
  ~IVector();
  void resize(const int &m);
  void assign(const int &m,const int &x);
  IVector& operator=(const IVector &v);
  int n;
  int *t;
};

class RMatrix
{
public:
  RMatrix();
  RMatrix(const int &nn,const int &mm);
  RMatrix(const int &nn,const int &mm,const double &z);
  RMatrix(const RMatrix &M);
  ~RMatrix();
  void resize(const int &nn,const int &mm);
  void assign(const int &nn,const int &mm,const double &z);
  RMatrix& operator=(const RMatrix &M);
  int n,m;
  double **t;
};

double sum(const RVector &V);

std::ostream &operator<<(std::ostream &out,const RVector &v);
std::ostream &operator<<(std::ostream &out,const RMatrix &M);

#endif // REALCONTAINERS_H_INCLUDED
