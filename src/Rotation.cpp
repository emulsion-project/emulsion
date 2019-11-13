#include "Rotation.h"

Rotation::Rotation() {}

Rotation::Rotation(const int &n, const double &a, const double &b, const double &g) {
  N = n;
  NM = (N+1)*(N+1)-1;
  alpha = a;
  beta = b;
  gamma = g;
  wigner.assign(N,beta);
}

void Rotation::setParameters(const int &n, const double &a, const double &b, const double &g) {
  N = n;
  NM = (N+1)*(N+1)-1;
  alpha = a;
  beta = b;
  gamma = g;
  wigner.assign(N,beta);
}

void Rotation::calcRotationMatrices() {
  Rd = Ri = Matrix(NM,NM,0.);
  int row = 0, col;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      col = 0;
      for (int l=1 ; l<=N ; l++) {
        for (int k=-l ; k<=l ; k++) {
          if (n == l) {
            Rd.t[row][col] = exp(i_*(double(m)*alpha+double(k)*gamma))*wigner.dlmn(k,m,n);
            Ri.t[col][row] = conj(Rd.t[row][col]);
          }
          col++;
        }
      }
      row++;
    }
  }
}

Matrix Rotation::getRd() {
  Matrix R(2*NM,2*NM);
  for (int i=0 ; i<NM ; i++) for (int j=0 ; j<NM ; j++) R.t[i][j] = R.t[i+NM][j+NM] = Rd.t[i][j];
  return R;
}

Matrix Rotation::getRi() {
  Matrix R(2*NM,2*NM);
  for (int i=0 ; i<NM ; i++) for (int j=0 ; j<NM ; j++) R.t[i][j] = R.t[i+NM][j+NM] = Ri.t[i][j];
  return R;
}

Vector Rotation::rotateVectorD(const Vector &v) {
  Vector vv(2*NM,0.);
  int row = 0;
  Complex D;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      for (int k=-n ; k<=n ; k++) {
        D = exp(i_*(double(m)*alpha+double(k)*gamma))*wigner.dlmn(k,m,n);
        vv.t[row] += D*v.t[ii(n,k)];
        vv.t[row+NM] += D*v.t[ii(n,k)+NM];
      }
      row++;
    }
  }
  return vv;
}

Vector Rotation::rotateVectorI(const Vector &v) {
  Vector vv(2*NM,0.);
  int row = 0;
  Complex D;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      for (int k=-n ; k<=n ; k++) {
        D = exp(-i_*(double(k)*alpha+double(m)*gamma))*wigner.dlmn(m,k,n);
        vv.t[row] += D*v.t[ii(n,k)];
        vv.t[row+NM] += D*v.t[ii(n,k)+NM];
      }
      row++;
    }
  }
  return vv;
}
