#include "Matrix.h"
#include <iostream>
using namespace std;

/*********************************
 *     Class constructors        *
 *********************************/

Matrix::Matrix() : n(0), m(0), t(NULL) {}

Matrix::Matrix(const int &nn, const int &mm) : n(nn), m(mm), t(nn>0 ? new Complex*[nn] : NULL) {
  int mn=mm*nn;
  if (t) t[0] = mn>0 ? new Complex[mn] : NULL;
  for (int i=1 ; i<nn ; i++) t[i] = t[i-1] + mm;
}

Matrix::Matrix(const int &nn, const int &mm, const Complex &x) : n(nn), m(mm), t(n>0 ? new Complex*[nn] : NULL) {
  int mn = mm*nn;
  if (t) t[0] = mn>0 ? new Complex[mn] : NULL;
  for (int i=1 ; i<nn ; i++) t[i] = t[i-1] + mm;
  for (int i=0 ; i<nn ; i++) for (int j=0 ; j<mm ; j++) t[i][j] = x;
}

Matrix::Matrix(const Matrix &M) : n(M.n), m(M.m), t(n>0 ? new Complex*[n] : NULL) {
  int mn = n*m;
  if (t != NULL) t[0] = mn>0 ? new Complex[mn] : NULL;
  for (int i=1 ; i<n ; i++) t[i] = t[i-1] + m;
  for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = M.t[i][j];
}

/*********************************
 *     Class destructor          *
 *********************************/

Matrix::~Matrix() {
  if (t != NULL) {
    delete[] (t[0]);
    delete[] (t);
  }
}

/*********************************
 *     Class member functions    *
 *********************************/

void Matrix::resize(const int &nn,const int &mm) {
  int mn;
  if (nn != n || mm != m) {
    if (t != NULL) {
      delete[] (t[0]);
      delete[] (t);
    }
    n = nn;
    m = mm;
    mn = n*m;
    t = n>0 ? new Complex*[n] : NULL;
    if (t) t[0] = mn>0 ? new Complex[mn] : NULL;
    for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
  }
}

void Matrix::assign(const int &nn,const int &mm,const Complex &x) {
  int mn;
  if (nn != n || mm != m) {
    if(t!=NULL) {
      delete[] (t[0]);
      delete[] (t);
    }
    n = nn;
    m = mm;
    mn = n*m;
    t = n>0 ? new Complex*[n] : NULL;
    if (t) t[0] = mn>0 ? new Complex[mn] : NULL;
    for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
  }
  for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = x;
}

Vector Matrix::col(const int &i) {
  Vector V(n);
  for(int j=0 ; j<n ; j++) V.t[j] = t[j][i];
  return V;
}

Vector Matrix::row(const int &i) {
  Vector V(m);
  for(int j=0 ; j<m ; j++) V.t[j] = t[i][j];
  return V;
}

/*********************************
 *     Class member operators    *
 *********************************/

Complex& Matrix::operator()(const int &i,const int &j) const { return t[i][j]; }
Complex& Matrix::operator()(const int &i,const int &j) { return t[i][j]; }

Matrix& Matrix::operator=(const Matrix &M) {
  if (this != &M) {
    int mn;
    if ((this->n != M.n) || (this->m != M.m)) {
      if (t != NULL) {
        delete[] (t[0]);
        delete[] (t);
      }
      n = M.n;
      m = M.m;
      t = n>0 ? new Complex*[n] : NULL;
      mn = n*m;
      if (t) t[0] = mn>0 ? new Complex[mn] : NULL;
      for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
    }
    for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = M.t[i][j];
  }
  return *this;
}

/*********************************
 *     Associated operators      *
 *********************************/

Matrix operator+(const Matrix &M,const Matrix &N) {
  Matrix X(M.n,M.m,0.);
  for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<M.m ; j++) X.t[i][j] = M.t[i][j] + N.t[i][j];
  return X;
}

Matrix operator-(const Matrix &M,const Matrix &N) {
  Matrix X(M.n,M.m);
  for (int i=0 ; i<M.n ; i++)  for (int j=0 ; j<M.m ; j++) X.t[i][j] = M.t[i][j] - N.t[i][j];
  return X;
}

Matrix operator*(const Matrix &M,const Matrix &N) {
  Matrix X(M.n,N.m,0.0);
  for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<N.m ; j++) for (int k=0 ; k<M.m ; k++) X.t[i][j] += M.t[i][k] * N.t[k][j];
  return X;
}

Matrix operator*(const Complex &z,const Matrix &M) {
  Matrix X = M;
  for (int i=0 ; i<X.n ; i++) for (int j=0 ; j<X.m ; j++) X.t[i][j] *= z;
  return X;
}

Matrix operator*(const Matrix &M,const Complex &z) {
  Matrix X = M;
  for (int i=0 ; i<X.n ; i++) for (int j=0 ; j<X.m ; j++) X.t[i][j] *= z;
  return X;
}

Vector operator*(const Matrix &M,const Vector &V) {
  Vector X(V.n,0.);
  if (M.m == V.n) {
    for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<V.n ; j++) X.t[i] += M.t[i][j] * V.t[j];
    return X;
  }
  return X;
}

Vector operator*(const Vector &V,const Matrix &M) {
  Vector X(V.n,0.);
  if (M.n == V.n) {
    for (int i=0 ; i<M.m ; i++) for (int j=0 ; j<V.n ; j++) X.t[i] += M.t[j][i] * V.t[j];
    return X;
  }
  return X;
}

std::ostream &operator<<(std::ostream &out,const Matrix &M) {
  for (int i=0 ; i<M.n ; i++) {
    for (int j=0 ; j<M.m ; j++) out << M.t[i][j] << "\t";
    out << "\n";
  }
  return out;
}

/*********************************
 *     Associated functions      *
 *********************************/

Vector mult(void *A,const Vector &V) {
  Vector X(V.n,0.0);
  Matrix M = *(Matrix*)A;
  if (M.m == V.n) {
    for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<V.n ; j++) X.t[i] += M.t[i][j] * V.t[j];
    return X;
  }
  else return X;
}

Complex min(const Matrix &M) {
  Complex x = M.t[0][0];
  for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<M.m ; j++) if (x > M.t[i][j]) x = M.t[i][j];
  return x;
}

Complex max(const Matrix &M) {
  Complex x = M.t[0][0];
  for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<M.m ; j++) if (x < M.t[i][j]) x = M.t[i][j];
  return x;
}

bool gaussj(Matrix &a,Matrix &b) {
  int i,icol=0,irow=0,j,k,l,ll,n=a.n,m=b.m;
  double big;
  Complex dum,pivinv;
  int *indxc,*indxr,*ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
  for (i=0 ; i<n ; i++) ipiv[i] = 0;
  for (i=0 ; i<n ; i++) {
    big = 0.;
    for (j=0 ; j<n ; j++) {
      if (ipiv[j] != 1) {
        for (k=0 ; k<n ; k++) {
          if (ipiv[k] == 0) {
            if (abs(a.t[j][k]) >= big) {
              big = abs(a.t[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0 ; l<n ; l++) swap(a.t[irow][l],a.t[icol][l]);
      for (l=0 ; l<m ; l++) swap(b.t[irow][l],b.t[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(a.t[icol][icol]==0.) return false;
    pivinv = 1./a.t[icol][icol];
    a.t[icol][icol] = 1.;
    for(l=0 ; l<n ; l++) a.t[icol][l] *= pivinv;
    for(l=0 ; l<m ; l++) b.t[icol][l] *= pivinv;
    for(ll=0 ; ll<n ; ll++) {
      if(ll != icol) {
          dum = a.t[ll][icol];
          a.t[ll][icol] = 0.0;
          for (l=0 ; l<n ; l++) a.t[ll][l] -= a.t[icol][l]*dum;
          for (l=0 ; l<m ; l++) b.t[ll][l] -= b.t[icol][l]*dum;
      }
    }
  }
  for (l=n-1 ; l>=0 ; l--) if (indxr[l] != indxc[l]) for (k=0 ; k<n ; k++) swap(a.t[k][indxr[l]],a.t[k][indxc[l]]);
  return true;
}

bool gaussj(Matrix &a,Vector &b) {
  int i,icol=0,irow=0,j,k,l,ll,n=a.n;
  icol = 0;
  double big;
  Complex dum,pivinv;
  int *indxc,*indxr,*ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
  for (i=0 ; i<n ; i++) ipiv[i] = 0;
  for (i=0 ; i<n ; i++) {
    big = 0.;
    for (j=0 ; j<n ; j++) {
      if (ipiv[j] != 1) {
        for (k=0 ; k<n ; k++) {
          if (ipiv[k] == 0) {
            if (abs(a.t[j][k]) >= big) {
              big = abs(a.t[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0 ; l<n ; l++) swap(a.t[irow][l],a.t[icol][l]);
      swap(b.t[irow],b.t[icol]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a.t[icol][icol] == 0.0) return false;
    pivinv = 1./a.t[icol][icol];
    a.t[icol][icol] = 1.0;
    for(l=0 ; l<n ; l++) a.t[icol][l] *= pivinv;
    b.t[icol] *= pivinv;
    for(ll=0 ; ll<n ; ll++) {
      if(ll != icol) {
        dum = a.t[ll][icol];
        a.t[ll][icol] = 0.0;
        for (l=0 ; l<n ; l++) a.t[ll][l] -= a.t[icol][l]*dum;
        b.t[ll] -= b.t[icol]*dum;
      }
    }
  }
  for (l=n-1 ; l>=0 ; l--) if (indxr[l] != indxc[l]) for (k=0 ; k<n ; k++) swap(a.t[k][indxr[l]],a.t[k][indxc[l]]);
  return true;
}

bool invGaussj(Matrix &a) {
  int i,icol=0,irow=0,j,k,l,ll,n=a.n;
  double big;
  Complex dum,pivinv;
  int *indxc,*indxr,*ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
  for (i=0 ; i<n ; i++) ipiv[i] = 0;
  for (i=0 ; i<n ; i++) {
    big = 0.;
    for (j=0 ; j<n ; j++) {
      if (ipiv[j] != 1) {
        for (k=0 ; k<n ; k++) {
          if (ipiv[k] == 0) {
            if (abs(a.t[j][k]) >= big) {
              big = abs(a.t[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    if (irow != icol)  for (l=0 ; l<n ; l++) swap(a.t[irow][l],a.t[icol][l]);
    indxr[i]=irow;
    indxc[i]=icol;
    if (a.t[icol][icol] == 0.) return false;
    pivinv = 1./a.t[icol][icol];
    a.t[icol][icol] = 1.;
    for(l=0 ; l<n ; l++) a.t[icol][l] *= pivinv;
    for(ll=0 ; ll<n ; ll++) {
      if(ll != icol) {
        dum = a.t[ll][icol];
        a.t[ll][icol] = 0.;
        for (l=0 ; l<n ; l++) a.t[ll][l] -= a.t[icol][l]*dum;
      }
    }
  }
  for (l=n-1 ; l>=0 ; l--) if (indxr[l] != indxc[l]) for (k=0 ; k<n ; k++) swap(a.t[k][indxr[l]],a.t[k][indxc[l]]);
  return true;
}

void hessenberg(Matrix &A) {
  int n = A.m, i, j, k, piv=0;
  Complex max = 0., ri;
  for (i=0 ; i<n-2 ; i++) {
    max = 0.;
    for (j=i+1 ; j<n ; j++) if (abs(max) < abs(A.t[j][i])) { max = A.t[j][i]; piv = j; }
    if (piv != i+1) {
      for (j=i ; j<n ; j++) swap(A.t[piv][j],A.t[i+1][j]);
      for (j=0 ; j<n ; j++) swap(A.t[j][piv],A.t[j][i+1]);
    }
    if (abs(max) > 0.) {
      for (j=i+2 ; j<n ; j++) {
        if (abs(A.t[j][i]) > 0.) {
          ri = A.t[j][i]/max;
          for (k=i ; k<n ; k++) A.t[j][k] -= ri*A.t[i+1][k];
          for (k=0 ; k<n ; k++) A.t[k][i+1] += ri*A.t[k][j];
        }
      }
    }
  }
}

void qr(Matrix &A) {
  hessenberg(A);
  int m = A.m, i, j, k, imin, imax;
  double a;
  Complex a1, a2, a3, a4, lambda1, lambda, ri;
  Vector c(m-1), s(m-1);
  imax = m-1;
  imin = 0;
  for (k=0 ; k<10*m ; k++) {
    while ((imax > imin) && (abs(A.t[imax][imax-1]) < 1e-15)) imax--;
    if (imax == imin) break;
    a1 = A.t[imax-1][imax-1];
    a2 = A.t[imax-1][imax];
    a3 = A.t[imax][imax-1];
    a4 = A.t[imax][imax];
    lambda1 = (a1+a4)*0.5 + sqrt((a1-a4)*(a1-a4)*0.25+a2*a3);
    lambda = (a1+a4)*0.5 - sqrt((a1-a4)*(a1-a4)*0.25+a2*a3);
    if (abs(a4-lambda1) < abs(a4-lambda)) lambda = lambda1;
    for (i=0 ; i<=imax ; i++) A.t[i][i] -= lambda;
    for (i=0 ; i<imax ; i++) {
      a1 = A.t[i][i];
      a2 = A.t[i+1][i];
      a = sqrt(norm(a1)+norm(a2));
      c.t[i] = a1/a;
      s.t[i] = a2/a;
      for (j=0 ; j<m ; j++) {
        a1 = A.t[i][j];
        a2 = A.t[i+1][j];
        A.t[i][j] = conj(c.t[i])*a1 + conj(s.t[i])*a2;
        A.t[i+1][j] = c.t[i]*a2 - s.t[i]*a1;
      }
    }
    for (i=0 ; i<imax ; i++) {
      for (j=0 ; j<m ; j++) {
        a1 = A.t[j][i];
        a2 = A.t[j][i+1];
        A.t[j][i] = c.t[i]*a1 + s.t[i]*a2;
        A.t[j][i+1] = conj(c.t[i])*a2 - conj(s.t[i])*a1;
      }
    }
    for (i=0 ; i<=imax ; i++) A.t[i][i] += lambda;
  }
  for (i=0 ; i<m-1 ; i++) {
    for (j=i ; j<m-1 ; j++) {
      ri = A.t[i][j+1]/(A.t[i][i] - A.t[j+1][j+1]);
      for (k=0 ; k<m ; k++) A.t[i][k] += ri*A.t[j+1][k];
      for (k=0 ; k<m ; k++) A.t[k][j+1] -= ri*A.t[k][i];
    }
  }
}

Vector roots(const Vector &AA) {
  int i, n = AA.n;
  Matrix A(n,n,0.);
  Vector R(n);
  for (i=0 ; i<n-1 ; i++) A.t[i+1][i] = 1.;
  for (i=0 ; i<n ; i++) A.t[i][n-1] = -AA.t[i];
  qr(A);
  for (i=0 ; i<n ; i++) R.t[i] = A.t[i][i];
  return R;
}

Vector poles(const Vector &Y,const Vector &X,const int &m) {
  int i, j, k, n = Y.n-2*m+1;
  Vector XX = X;
  Complex mil = abs(XX.t[XX.n-1] - XX.t[0]);
  for (i=0 ; i<XX.n ; i++) XX.t[i] -= mil;
  Matrix A(m,m,0.);
  Vector V(m), B(m);
  for (i=0 ; i<m+n ; i++) {
    for (j=0 ; j<m ; j++) V.t[j] = Y.t[i+j];
    for (k=0 ; k<m+n ; k++) if (k != i) for (j=0 ; j<m ; j++) V.t[j] /= XX.t[i+j] - XX.t[k+j];
    for (j=0 ; j<m ; j++) {
      for (k=0 ; k<m ; k++) {
        A.t[j][k] += V.t[j];
        V.t[j] *= XX.t[i+j];
      }
      B.t[j] -= V.t[j];
    }
  }
  gaussj(A,B);
  B = roots(B);
  for (i=0 ; i<B.n ; i++) B.t[i] += mil;
  return B;
}

Matrix polesAmp(const Vector &Y,const Vector &X,const int &m) {
  int i, j, k, n = Y.n-2*m+1;
  Vector XX = X;
  //Complex a;
  Complex mil = abs(XX.t[XX.n-1] - XX.t[0]),a;
  for (i=0 ; i<XX.n ; i++) XX.t[i] -= mil;
  Matrix A(m,m,0.), PA(m,2,0.);
  Vector V(m), B(m);
  for (i=0 ; i<m+n ; i++) {
    for (j=0 ; j<m ; j++) V.t[j] = Y.t[i+j];
    for (k=0 ; k<m+n ; k++) if (k != i) for (j=0 ; j<m ; j++) V.t[j] /= XX.t[i+j] - XX.t[k+j];
    for (j=0 ; j<m ; j++) {
      for (k=0 ; k<m ; k++) {
        A.t[j][k] += V.t[j];
        V.t[j] *= XX.t[i+j];
      }
      B.t[j] -= V.t[j];
    }
  }
  gaussj(A,B);
  B = roots(B);


  for (int i=0 ; i<m ; i++) {
    PA.t[i][1] = 0.;
    for (int j=0 ; j<Y.n ; j++) {
      a = Y.t[j];
      for (int ii=0 ; ii<m ; ii++) a *= (XX.t[j] - B.t[ii]);
      for (int ii=0 ; ii<Y.n ; ii++) if (ii != j) a *= (B.t[i] - XX.t[ii])/(XX.t[j] - XX.t[ii]);
      PA.t[i][1] += a;
    }
    for (int j=0 ; j<m ; j++) if (j != i) PA.t[i][1] /= B.t[i] - B.t[j];
    //PA.t[i][1] *= B.t[i];
  }

  for (i=0 ; i<B.n ; i++) PA.t[i][0] = B.t[i]+mil;

  return PA;
}

/*Matrix polesAmp(const Vector &Y,const Vector &X,const int &m) {
  int i, j, k, N = Y.n/m;
  Vector XX = X;
  //Complex a;
  Complex mil = abs(XX.t[XX.n-1] - XX.t[0]),a;
  for (i=0 ; i<XX.n ; i++) XX.t[i] -= mil;
  Matrix A(m,m,0.), PA(m,2,0.);
  Vector V(m), B(m);
  for (i=0 ; i<m ; i++) {
    for (j=0 ; j<m ; j++) {
      for (int ii=0 ; ii<N ; ii++) {
        a = pow(XX.t[ii+N*i],j)*Y.t[ii+N*i];
        for (int l=0 ; l<N ; l++) {
          if (l != ii) a /= (XX.t[ii+N*i]-XX.t[l+N*i]);
        }
        A.t[i][j] += a;
      }
    }
    for (int ii=0 ; ii<N ; ii++) {
        a = pow(XX.t[ii+N*i],m)*Y.t[ii+N*i];
        for (int l=0 ; l<N ; l++) {
          if (l != ii) a /= (XX.t[ii+N*i]-XX.t[l+N*i]);
        }
        B.t[i] -= a;
      }
  }
  gaussj(A,B);
  B = roots(B);

  for (int i=0 ; i<m ; i++) {
    PA.t[i][1] = 0.;
    for (int j=0 ; j<Y.n ; j++) {
      a = Y.t[j];
      for (int ii=0 ; ii<m ; ii++) a *= (XX.t[j] - B.t[ii]);
      for (int ii=0 ; ii<Y.n ; ii++) if (ii != j) a *= (B.t[i] - XX.t[ii])/(XX.t[j] - XX.t[ii]);
      PA.t[i][1] += a;
    }
    for (int j=0 ; j<m ; j++) if (j != i) PA.t[i][1] /= B.t[i] - B.t[j];
    //PA.t[i][1] *= B.t[i];
  }

  for (i=0 ; i<B.n ; i++) PA.t[i][0] = B.t[i]+mil;

  return PA;
}*/

Vector zeros(const Vector &Y,const Vector &X,const int &m) {
  int i, j, k, n = Y.n-2*m+1;
  Vector XX = X;
  Complex mil = XX.t[XX.n-1] - XX.t[0];
  for (i=0 ; i<XX.n ; i++) XX.t[i] -= mil;
  Matrix A(m,m,0.);
  Vector V(m), B(m);
  for (i=0 ; i<m+n ; i++) {
    for (j=0 ; j<m ; j++) V.t[j] = 1./Y.t[i+j];
    for (k=0 ; k<m+n ; k++) if (k != i) for (j=0 ; j<m ; j++) V.t[j] /= XX.t[i+j] - XX.t[k+j];
    for (j=0 ; j<m ; j++) {
      for (k=0 ; k<m ; k++) {
        A.t[j][k] += V.t[j];
        V.t[j] *= XX.t[i+j];
      }
      B.t[j] -= V.t[j];
    }
  }
  gaussj(A,B);
  B = roots(B);
  for (i=0 ; i<B.n ; i++) B.t[i] += mil;
  return B;
}

void gmres(Vector &X, void *Struct, Vector (*mult)(void *, const Vector &), const double &tol) {
  const int kmax = 1000;
  int i, k, it, itt=0, n=X.n;
  double res;
  Complex *h1,*h2,aux,ss[kmax],cc[kmax],gg[kmax+1],y[kmax];
  Vector XX, B, *V[kmax+1], *H[kmax+1];
  XX = B = X;
  do {
    it = 0;
    for (i=0; i<n ; i++) X.t[i] = 0.;
    for (i=0; i<kmax+1 ; i++) gg[i] = 0.;
    V[0] = new Vector(n);
    gg[0] = res = norm(XX);
    if (res > tol) res = 1./res;
    for (i=0 ; i<n ; i++) V[0]->t[i] = res*XX.t[i];
    do {
      V[it+1] = new Vector(n);
      H[it] = new Vector(it+2);
      *V[it+1] = mult(Struct,*V[it]);
      for (k=0 ; k<it+1 ; k++) {
        H[it]->t[k] = aux = (*V[k])*(*V[it+1]);
        for (i=0 ; i<n ; i++) V[it+1]->t[i] -= aux*V[k]->t[i];
      }
      H[it]->t[it+1] = res = norm(*V[it+1]);
      if (fabs(res)>1e-10) {
        res = 1./res;
        for (i=0 ; i<n ; i++) V[it+1]->t[i] *= res;
      }
      h1 = h2 = H[it]->t;
      h2++;
      for (k=0 ; k<it ; k++, h1++, h2++) {
        *h1 = cc[k]*(aux = *h1) + ss[k]*(*h2);
        *h2 = -conj(ss[k])*aux + cc[k]*(*h2);
      }
      h1 = h2 = H[it]->t + it;
      h2++;
      if (abs(*h2)<1e-10) {
        ss[it] = 0.;
        cc[it] = 1.;
      }
      else if (abs(*h1)<1e-10) {
        ss[it] = conj(*h2)/abs(*h2);
        cc[it] = 0.;
      }
      else {
        res = 1./sqrt(norm(*h1)+norm(*h2));
        cc[it] = res*abs(*h1);
        ss[it] = res*conj(*h2)*(*h1)/abs(*h1);
      }
      *h1 = cc[it]*(aux = *h1) + ss[it]*(*h2);
      *h2 = -conj(ss[it])*aux + cc[it]*(*h2);
      gg[it+1] = -gg[it]*conj(ss[it]);
      gg[it] *= cc[it];
      res = abs(gg[it+1]);
      cout << it << " residual " << res << endl;
      it++;
      itt++;
    } while ((it<kmax) && (res > tol));
    cout << itt << " residual: " << res << endl;
    for (k=it-1 ; k>=0 ; k--) {
      y[k] = gg[k];
      for (i=it-1 ; i>k ; i--) y[k] -= y[i]*H[i]->t[k];
      y[k] /= H[k]->t[k];
    }
    for (k=0 ; k<it ; k++) for (i=0 ; i<n ; i++) X.t[i] += y[k]*V[k]->t[i];
    for (k=0 ; k<it ; k++) {
      delete V[k];
      delete H[k];
    }
    delete V[it];
    XX = B - mult(Struct,X);
  } while (res > tol);
}

Vector gmres(const Vector &B, const Matrix &A, const double &tol) {
  const int kmax = 1000;
  int i, k, it=0, n=B.n;
  double res;
  Complex *h1,*h2,aux,ss[kmax],cc[kmax],gg[kmax+1],y[kmax];
  Vector X(B), *V[kmax+1], *H[kmax+1];
  for (i=0; i<kmax+1 ; i++) gg[i] = 0.;
  V[0] = new Vector(n);
  gg[0] = res = norm(X);
  if (res > tol) res = 1./res;
  for (i=0 ; i<n ; i++) V[0]->t[i] = res*X.t[i];
  do {
    V[it+1] = new Vector(n);
    H[it] = new Vector(it+2);
    *V[it+1] = A*(*V[it]);
    for (k=0 ; k<it+1 ; k++) {
      H[it]->t[k] = aux = (*V[k])*(*V[it+1]);
      for (i=0 ; i<n ; i++) V[it+1]->t[i] -= aux*V[k]->t[i];
    }
    H[it]->t[it+1] = res = norm(*V[it+1]);
    if (fabs(res)>tol) {
      res = 1./res;
      for (i=0 ; i<n ; i++) V[it+1]->t[i] *= res;
    }
    h1 = h2 = H[it]->t;
    h2++;
    for (k=0 ; k<it ; k++, h1++, h2++) {
      *h1 = cc[k]*(aux = *h1) + ss[k]*(*h2);
      *h2 = -conj(ss[k])*aux + cc[k]*(*h2);
    }
    h1 = h2 = H[it]->t + it;
    h2++;
    if (abs(*h2)<1e-10) {
      ss[it] = 0.;
      cc[it] = 1.;
    }
    else if (abs(*h1)<1e-10) {
      ss[it] = conj(*h2)/abs(*h2);
      cc[it] = 0.;
    }
    else {
      res = 1./sqrt(norm(*h1)+norm(*h2));
      cc[it] = res*abs(*h1);
      ss[it] = res*conj(*h2)*(*h1)/abs(*h1);
    }
    *h1 = cc[it]*(aux = *h1) + ss[it]*(*h2);
    *h2 = -conj(ss[it])*aux + cc[it]*(*h2);
    gg[it+1] = -gg[it]*conj(ss[it]);
    gg[it] *= cc[it];
    res = abs(gg[it+1]);
    //cout << it << " residual " << res << endl;
    it++;
  } while ((it<kmax) && (res > tol));
  for (k=it-1 ; k>=0 ; k--) {
    y[k] = gg[k];
    for (i=it-1 ; i>k ; i--) y[k] -= y[i]*H[i]->t[k];
    y[k] /= H[k]->t[k];
  }
  for (i=0; i<n ; i++) X.t[i] = 0.;
  for (k=0 ; k<it ; k++) for (i=0 ; i<n ; i++) X.t[i] += y[k]*V[k]->t[i];
  for (k=0 ; k<it ; k++) {
    delete V[k];
    delete H[k];
  }
  delete V[it];
  return X;
}

Matrix gmres(const Matrix &B, const Matrix &A, const double &tol) {
  Matrix X(B.n,B.m,0.);
  Vector BB(B.n),XX;
  for (int i=0 ; i<B.m ; i++) {
    for (int j=0 ; j<B.n ; j++) BB.t[j] = B.t[j][i];
    XX = gmres(BB,A,tol);
    for (int j=0 ; j<B.n ; j++) X.t[j][i] = XX.t[j];
  }
  return X;
}

Vector polyFit(const Matrix &XY, const int &n) {
  int N = XY.n;
  Matrix M(n+1,n+1,0.);
  Vector V(n+1,0.),X(n+1,0.);

  for (int i=0 ; i<=n ; i++) {
      cout << i << endl;
    for (int j=0 ; j<=n ; j++) {
      for (int nn=0 ; nn<N ; nn++) {
        M.t[i][j] += pow(XY.t[nn][0],double(i+j));
      }
    }
    for (int nn=0 ; nn<N ; nn++) {
      V.t[i] += XY.t[nn][1]*pow(XY.t[nn][0],double(i));
    }
  }
  invGaussj(M);
  X = M*V;
  return X;
}

Matrix transpose(const Matrix &M) {
  Matrix A(M.n,M.m);
  for (int i=0 ; i<M.n ; i++) for (int j=0 ; j<M.m ; j++) A.t[i][j] = M.t[j][i];
  return A;
}
