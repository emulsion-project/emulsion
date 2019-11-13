#include "Translation.h"
#include <iostream>
using namespace std;

Translation::Translation() {}

Translation::Translation(const int &typ, const bool &dir, const int &n, const Complex kk, const double &xx, const double &yy, const double &zz) {
  setParameters(typ,dir,n,kk,xx,yy,zz);
}

void Translation::setParameters(const int &typ, const bool &dir, const int &n, const Complex kk, const double &xx, const double &yy, const double &zz) {
  type = typ; direct = dir;
  N = n; NM = (N+1)*(N+1)-1;
  x = xx; y = yy; z = zz; k = kk;
  r = sqrt(x*x+y*y+z*z); t = acos(z/r); p = atan2(y,x);
  if (r==0.) t = 0.;
  if (y==0. && x==0.) p = 0.;
  cg.assign(N);
  bes.assign(2*N,k*r);
  if (direct) calcMatrix();
  else calcMatrixZ();
}

void Translation::calcMatrix() {
  Ad = Ai = Bd = Bi = Matrix(NM,NM);
  Complex A,B,mult,sb;
  leg.assign(2*N,cos(t));
  int col,row;
  row = 0;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      col = 0;
      for (int l=1 ; l<=N ; l++) {
        for (int k=-l ; k<=l ; k++) {
          A = B = 0.;
          for (int w=abs(n-l) ; w<=n+l ; w++) {
            sb = type == 1 ? bes.jn(w) : bes.h1n(w);
            sb *= powi(w)/sqrt(double(w+w+1))*leg.Pmn(w,m-k)*cg.coef(w,n,m,l,-k)*cg.coef(w,n,1,l,-1);
            if ((n+l-w)%2 == 0) A += sb;
            else B += sb;
          }
          mult = powm1(k+1)*powi(l-n)*sqrt(2.*double(n+n+1)*double(l+l+1))*exp(i_*double(m-k)*p);
          A *= mult; B *= mult;
          Ad.t[row][col] = A; Bd.t[row][col] = B;
          Ai.t[row][col] = powm1(n+l)*A; Bi.t[row][col] = powm1(n+l+1)*B;
          col++;
        }
      }
      row++;
    }
  }
}

void Translation::calcMatrixZ() {
  Ad = Ai = Bd = Bi = Matrix(NM,NM);
  Complex A,B,mult,sb;
  wigner.assign(N,t);
  int col,row;
  row = 0;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      col = 0;
      for (int l=1 ; l<=N ; l++) {
        for (int k=-l ; k<=l ; k++) {
          if (m == k) {
            A = B = 0.;
            for (int w=abs(n-l) ; w<=n+l ; w++) {
              sb = (type == 1) ? bes.jn(w) : bes.h1n(w);
              sb *= powi(w)*cg.coef(w,n,m,l,-m)*cg.coef(w,n,1,l,-1);
              if ((n+l-w)%2 == 0) A += sb;
              else B += sb;
            }
            mult = powm1(m+1)*powi(l-n)*sqrt(double(n+n+1)*double(l+l+1));
            A *= mult; B *= mult;
            Ad.t[row][col] = A; Bd.t[row][col] = B;
            Ai.t[row][col] = powm1(n+l)*A; Bi.t[row][col] = powm1(n+l+1)*B;
          }
          else Ad.t[row][col] = Bd.t[row][col] = Ai.t[row][col] = Bi.t[row][col] = 0.;
          col++;
        }
      }
      row++;
    }
  }
}

Vector Translation::translate(const Vector &v, const int &direction) {
  Vector vv(2*NM,0.);
  Vector w1(2*NM,0.);
  if (direct) {
    for (int i=0 ; i<NM ; i++) {
      for (int j=0 ; j<NM ; j++) {
        if (direction == 0) {
          vv.t[i] += v.t[j]*Ad.t[i][j] + v.t[j+NM]*Bd.t[i][j];
          vv.t[i+NM] += v.t[j]*Bd.t[i][j] + v.t[j+NM]*Ad.t[i][j];
        }
        else {
          vv.t[i] += v.t[j]*Ai.t[i][j] + v.t[j+NM]*Bi.t[i][j];
          vv.t[i+NM] += v.t[j]*Bi.t[i][j] + v.t[j+NM]*Ai.t[i][j];
        }
      }
    }
  }
  else {
    int row=0, col,n,l,m,k;
    Complex D;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        for (k=-n ; k<=n ; k++) {
          D = exp(-i_*double(k)*p)*wigner.dlmn(k,m,n);
          vv.t[row] += D*v.t[ii(n,k)]; vv.t[row+NM] += D*v.t[ii(n,k)+NM];
        }
        row++;
      }
    }
    row = 0;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        col = 0;
        for (l=1 ; l<=N ; l++) {
          for (k=-l ; k<=l ; k++) {
            if (k == m) {
              if (direction == 0) {
                w1.t[row] += vv.t[col]*Ad.t[row][col] + vv.t[col+NM]*Bd.t[row][col];
                w1.t[row+NM] += vv.t[col]*Bd.t[row][col] + vv.t[col+NM]*Ad.t[row][col];
              }
              else {
                w1.t[row] += vv.t[col]*Ai.t[row][col] + vv.t[col+NM]*Bi.t[row][col];
                w1.t[row+NM] += vv.t[col]*Bi.t[row][col] + vv.t[col+NM]*Ai.t[row][col];
              }
            }
            col++;
          }
        }
        row++;
      }
    }
    vv.assign(2*NM,0.);
    row = 0;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        for (k=-n ; k<=n ; k++) {
          D = exp(i_*(double(m)*p))*wigner.dlmn(m,k,n);
          vv.t[row] += D*w1.t[ii(n,k)]; vv.t[row+NM] += D*w1.t[ii(n,k)+NM];
        }
        row++;
      }
    }
  }
  return vv;
}

Vector Translation::translateT(const Vector &v, const int &direction) {
  Vector vv(2*NM,0.);
  if (direct) {
    for (int i=0 ; i<NM ; i++) {
      for (int j=0 ; j<NM ; j++) {
        if (direction == 0) {
          vv.t[i] += v.t[j]*Ad.t[j][i] + v.t[j+NM]*Bd.t[j][i];
          vv.t[i+NM] += v.t[j]*Bd.t[j][i] + v.t[j+NM]*Ad.t[j][i];
        }
        else {
          vv.t[i] += v.t[j]*Ai.t[j][i] + v.t[j+NM]*Bi.t[j][i];
          vv.t[i+NM] += v.t[j]*Bi.t[j][i] + v.t[j+NM]*Ai.t[j][i];
        }
      }
    }
  }
  else {
    int row=0, col,n,l,m,k;
    Vector w1(2*NM,0.);
    Complex D;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        for (k=-n ; k<=n ; k++) {
          D = exp(i_*(double(k)*p))*wigner.dlmn(k,m,n);
          vv.t[row] += D*v.t[ii(n,k)]; vv.t[row+NM] += D*v.t[ii(n,k)+NM];
        }
        row++;
      }
    }
    col = 0;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        row = 0;
        for (l=1 ; l<=N ; l++) {
          for (k=-l ; k<=l ; k++) {
            if (k == m) {
              if (direction == 0) {
                w1.t[row] += vv.t[col]*Ad.t[col][row] + vv.t[col+NM]*Bd.t[col][row];
                w1.t[row+NM] += vv.t[col]*Bd.t[col][row] + vv.t[col+NM]*Ad.t[col][row];
              }
              else {
                w1.t[row] += vv.t[col]*Ai.t[col][row] + vv.t[col+NM]*Bi.t[col][row];
                w1.t[row+NM] += vv.t[col]*Bi.t[col][row] + vv.t[col+NM]*Ai.t[col][row];
              }
            }
            row++;
          }
        }
        col++;
      }
    }
    vv.assign(2*NM,0.);
    row = 0;
    for (n=1 ; n<=N ; n++) {
      for (m=-n ; m<=n ; m++) {
        for (k=-n ; k<=n ; k++) {
          D = exp(-i_*(double(m)*p))*wigner.dlmn(m,k,n);
          vv.t[row] += D*w1.t[ii(n,k)]; vv.t[row+NM] += D*w1.t[ii(n,k)+NM];
        }
        row++;
      }
    }
  }
  return vv;
}

Matrix Translation::getMatrix(const int &direction) {
  Matrix T(2*NM,2*NM);
  for (int i=0 ; i<NM ; i++) {
    for (int j=0 ; j<NM ; j++) {
      if (direction == 0) {
        T.t[i][j] = T.t[i+NM][j+NM] = Ad.t[i][j];
        T.t[i+NM][j] = T.t[i][j+NM] = Bd.t[i][j];
      }
      else {
        T.t[i][j] = T.t[i+NM][j+NM] = Ai.t[i][j];
        T.t[i+NM][j] = T.t[i][j+NM] = Bi.t[i][j];
      }
    }
  }
  return T;
}

Matrix Translation::getMatrixT(const int &direction) {
  Matrix T(2*NM,2*NM);
  for (int i=0 ; i<NM ; i++) {
    for (int j=0 ; j<NM ; j++) {
      if (direction == 0) {
        T.t[j][i] = T.t[j+NM][i+NM] = Ad.t[i][j];
        T.t[j+NM][i] = T.t[j][i+NM] = Bd.t[i][j];
      }
      else {
        T.t[j][i] = T.t[j+NM][i+NM] = Ai.t[i][j];
        T.t[j+NM][i] = T.t[j][i+NM] = Bi.t[i][j];
      }
    }
  }
  return T;
}
