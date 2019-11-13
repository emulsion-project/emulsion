#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "Vector.h"

/*********************************
 *     Class declaration         *
 *********************************/

class Matrix
{
public:
  Matrix();
  Matrix(const int &nn,const int &mm);
  Matrix(const int &nn,const int &mm,const Complex &z);
  Matrix(const Matrix &M);
  ~Matrix();

  void resize(const int &nn,const int &mm);
  void assign(const int &nn,const int &mm,const Complex &z);
  Vector col(const int &i);
  Vector row(const int &i);

  Matrix& operator=(const Matrix &M);
  Complex& operator()(const int &i,const int &j);
  Complex& operator()(const int &i,const int &j) const;

  int n,m;
  Complex **t;
};

/*********************************
 *     Associated operators      *
 *********************************/

Matrix operator+(const Matrix &M,const Matrix &N);
Matrix operator-(const Matrix &M,const Matrix &N);
Matrix operator*(const Matrix &M,const Matrix &N);
Matrix operator*(const Complex &z,const Matrix &M);
Matrix operator*(const Matrix &M,const Complex &z);
Vector operator*(const Matrix &M,const Vector &V);
Vector operator*(const Vector &V,const Matrix &M);
std::ostream &operator<<(std::ostream &out,const Matrix &M);

/*********************************
 *     Associated functions      *
 *********************************/

Vector mult(void*,const Vector &V);
Complex min(const Matrix &M);
Complex max(const Matrix &M);
bool gaussj(Matrix &a, Matrix &b);
bool gaussj(Matrix &a,Vector &b);
bool invGaussj(Matrix &a);
void hessenberg(Matrix &A);
void qr(Matrix &A);
Vector roots(const Vector &AA);
Vector poles(const Vector &Y,const Vector &X,const int &m);
Matrix polesAmp(const Vector &Y,const Vector &X,const int &m);
Vector poleAmplitudes(const Vector &x, const Matrix &y);
Vector zeros(const Vector &Y,const Vector &X,const int &m);
void gmres(Vector &X, void *Struct, Vector (*mult)(void *, const Vector &), const double &tol);

Vector gmres(const Vector &B, const Matrix &A, const double &tol);
Matrix gmres(const Matrix &B, const Matrix &A, const double &tol);
Vector polyFit(const Matrix &XY, const int &n);
Matrix transpose(const Matrix &M);

#endif // MATRIX_H_INCLUDED
