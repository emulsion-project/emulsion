#include "RealContainers.h"

/*********************************
 *     Class constructors        *
 *********************************/

RVector::RVector() : n(0), t(NULL) {}
RVector::RVector(const int &n) : n(n), t(n>0? new double[n] : NULL) {}
RVector::RVector(const int &n,const double &x) : n(n), t(n>0 ? new double[n] : NULL) { for (int i=0 ; i<n ; i++) t[i] = x; }
RVector::RVector(const RVector &v) : n(v.n), t(n>0 ? new double[n] : NULL) { for (int i=0 ; i<n ; i++) t[i] = v.t[i]; }

IVector::IVector() : n(0), t(NULL) {}
IVector::IVector(const int &n) : n(n), t(n>0? new int[n] : NULL) {}
IVector::IVector(const int &n,const int &x) : n(n), t(n>0 ? new int[n] : NULL) { for (int i=0 ; i<n ; i++) t[i] = x; }
IVector::IVector(const IVector &v) : n(v.n), t(n>0 ? new int[n] : NULL) { for (int i=0 ; i<n ; i++) t[i] = v.t[i]; }

RMatrix::RMatrix() : n(0), m(0), t(NULL) {}

RMatrix::RMatrix(const int &nn, const int &mm) : n(nn), m(mm), t(nn>0 ? new double*[nn] : NULL) {
  int mn=mm*nn;
  if (t) t[0] = mn>0 ? new double[mn] : NULL;
  for (int i=1 ; i<nn ; i++) t[i] = t[i-1] + mm;
}

RMatrix::RMatrix(const int &nn, const int &mm, const double &x) : n(nn), m(mm), t(n>0 ? new double*[nn] : NULL) {
  int mn = mm*nn;
  if (t) t[0] = mn>0 ? new double[mn] : NULL;
  for (int i=1 ; i<nn ; i++) t[i] = t[i-1] + mm;
  for (int i=0 ; i<nn ; i++) for (int j=0 ; j<mm ; j++) t[i][j] = x;
}

RMatrix::RMatrix(const RMatrix &M) : n(M.n), m(M.m), t(n>0 ? new double*[n] : NULL) {
  int mn = n*m;
  if (t != NULL) t[0] = mn>0 ? new double[mn] : NULL;
  for (int i=1 ; i<n ; i++) t[i] = t[i-1] + m;
  for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = M.t[i][j];
}

/*********************************
 *     Class destructor          *
 *********************************/

RVector::~RVector() { if(t!=NULL) delete[] (t); }

IVector::~IVector() { if(t!=NULL) delete[] (t); }

RMatrix::~RMatrix() {
  if (t != NULL) {
    delete[] (t[0]);
    delete[] (t);
  }
}

/*********************************
 *     Class member operators    *
 *********************************/

RVector& RVector::operator=(const RVector &v) {
  if (this != &v) {
    if (n != v.n) {
      if (t !=	NULL) delete [] (t);
      n = v.n;
      t = n>0 ? new double[n] : NULL;
    }
    for (int i=0 ; i<n ; i++) t[i] = v.t[i];
  }
  return *this;
}

IVector& IVector::operator=(const IVector &v) {
  if (this != &v) {
    if (n != v.n) {
      if (t !=	NULL) delete [] (t);
      n = v.n;
      t = n>0 ? new int[n] : NULL;
    }
    for (int i=0 ; i<n ; i++) t[i] = v.t[i];
  }
  return *this;
}

RMatrix& RMatrix::operator=(const RMatrix &M) {
  if (this != &M) {
    int mn;
    if ((this->n != M.n) || (this->m != M.m)) {
      if (t != NULL) {
        delete[] (t[0]);
        delete[] (t);
      }
      n = M.n;
      m = M.m;
      t = n>0 ? new double*[n] : NULL;
      mn = n*m;
      if (t) t[0] = mn>0 ? new double[mn] : NULL;
      for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
    }
    for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = M.t[i][j];
  }
  return *this;
}


/*********************************
 *     Class member functions    *
 *********************************/

void RVector::resize(const int &m) {
  if (m != n) {
    if (t != NULL) delete[] (t);
    n = m;
    t = n>0 ? new double[n] : NULL;
  }
}

void RVector::assign(const int &m,const double &z) {
  if (m != n) {
    if (t != NULL) delete[] (t);
    n = m;
    t = n>0 ? new double[n] : NULL;
  }
  for (int i=0 ; i<n ; i++) t[i] = z;
}

void RMatrix::resize(const int &nn,const int &mm) {
  int mn;
  if (nn != n || mm != m) {
    if (t != NULL) {
      delete[] (t[0]);
      delete[] (t);
    }
    n = nn;
    m = mm;
    mn = n*m;
    t = n>0 ? new double*[n] : NULL;
    if (t) t[0] = mn>0 ? new double[mn] : NULL;
    for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
  }
}

void RMatrix::assign(const int &nn,const int &mm,const double &x) {
  int mn;
  if (nn != n || mm != m) {
    if(t!=NULL) {
      delete[] (t[0]);
      delete[] (t);
    }
    n = nn;
    m = mm;
    mn = n*m;
    t = n>0 ? new double*[n] : NULL;
    if (t) t[0] = mn>0 ? new double[mn] : NULL;
    for (int i=1 ; i<n; i++) t[i] = t[i-1] + m;
  }
  for (int i=0 ; i<n ; i++) for (int j=0 ; j<m ; j++) t[i][j] = x;
}


std::ostream &operator<<(std::ostream &out,const RVector &v) {
  for(int i=0 ; i<v.n ; i++) out << v.t[i] << "\n";
  return out;
}

std::ostream &operator<<(std::ostream &out,const RMatrix &M) {
  for (int i=0 ; i<M.n ; i++) {
    for (int j=0 ; j<M.m ; j++) out << M.t[i][j] << "\t";
    out << "\n";
  }
  return out;
}

double sum(const RVector &V) {
  double z=0.;
  for (int i=0 ; i<V.n ; i++) z += V.t[i];
  return z;
}
