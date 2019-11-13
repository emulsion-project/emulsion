#include "WignerFunctions.h"
#include <iostream>
using namespace std;

WignerFunctions::WignerFunctions() {
  N = 0;
  coefs = NULL;
}

WignerFunctions::~WignerFunctions() {
  if (coefs != NULL) {
    for (int i=0 ; i<N+N+1 ; i++) delete [] (coefs[i]);
    delete [] (coefs);
  }
}

WignerFunctions::WignerFunctions(const int &nn, const double &thet) {
  N = nn;
  theta = thet;
  x = cos(theta);
  int n,m,l,row,col;
  double d,dn,dnm1,ee;
  coefs = new RVector*[N+N+1];
  for (int r=0 ; r<N+N+1 ; r++) coefs[r] = new RVector[N+N+1];
  row = col = 0;
  double lns[N+N+1];
  lns[0] = lns[1] = 0.;
  for (int i=1 ; i<N+N+1 ; i++) lns[i+1] = lns[i] + log(double(i+1));

  for (m=-N ; m<=N ; m++) {
    col = 0;
    for (l=-N ; l<=N ; l++) {
      n0 = max(abs(m),abs(l));
      coefs[row][col] = RVector(N-n0+1);
      ee = l<m ? powm1(m-l) : 1.;
      coefs[row][col].t[0] = d = ee*pow(2.,-double(n0)) * sqrt(exp(lns[n0+n0]-lns[int(abs(m-l))]-lns[int(abs(m+l))])) * pow(1.-x,double(abs(m-l))*0.5)*pow(1.+x,double(abs(m+l))*0.5);
      dn = d;
      dnm1 = 0.;
      n = 1;
      for (n0 = n0 ; n0<N ; n0++) {
        if (n0==0) coefs[row][col].t[n] = d = double(n0+n0+1)*x*dn - dnm1;
        else d = 1./(double(n0) * sqrt(double((n0+1)*(n0+1))-double(m*m)) * sqrt(double((n0+1)*(n0+1))-double(l*l)))
          * (double(n0+n0+1)*(double(n0*n0+n0)*x-double(m*l))*dn - double(n0+1)*sqrt(double(n0*n0-m*m))*sqrt(double(n0*n0-l*l))*dnm1);
        dnm1 = dn;
        coefs[row][col].t[n] = dn = d;
        n++;
      }
      col++;
    }
    row++;
  }
}


void WignerFunctions::assign(const int &nn, const double &thet) {
  if (coefs != NULL) {
    for (int i=0 ; i<N+N+1 ; i++) delete [] (coefs[i]);
    delete [] (coefs);
  }
  N = nn;
  theta = thet;
  x = cos(theta);
  int n,m,l,row,col;
  double d,dn,dnm1,ee;
  coefs = new RVector*[N+N+1];
  for (int r=0 ; r<N+N+1 ; r++) coefs[r] = new RVector[N+N+1];
  row = col = 0;
  double lns[N+N+1];
  lns[0] = lns[1] = 0.;
  for (int i=1 ; i<N+N+1 ; i++) lns[i+1] = lns[i] + log(double(i+1));
  for (m=-N ; m<=N ; m++) {
    col = 0;
    for (l=-N ; l<=N ; l++) {
      n0 = max(abs(m),abs(l));
      coefs[row][col] = RVector(N-n0+1,0.);
      ee = l<m ? powm1(m-l) : 1.;
      coefs[row][col].t[0] = dn = ee*pow(2.,-double(n0)) * sqrt(exp(lns[n0+n0]-lns[int(abs(m-l))]-lns[int(abs(m+l))])) * pow(1.-x,double(abs(m-l))*0.5)*pow(1.+x,double(abs(m+l))*0.5);
      dnm1 = 0.;
      n = 1;
      for (n0 = n0 ; n0<N ; n0++) {
        if (n0==0) coefs[row][col].t[n] = d = double(n0+n0+1)*x*dn - dnm1;
        else coefs[row][col].t[n] = d = 1./(double(n0) * sqrt(double((n0+1)*(n0+1))-double(m*m)) * sqrt(double((n0+1)*(n0+1))-double(l*l)))
          * (double(n0+n0+1)*(double(n0*n0+n0)*x-double(m*l))*dn - double(n0+1)*sqrt(double(n0*n0-m*m))*sqrt(double(n0*n0-l*l))*dnm1);
        dnm1 = dn;
        dn = d;
        n++;
      }
      col++;
    }
    row++;
  }
}

double WignerFunctions::dlmn(const int &l, const int &m, const int &n) {
  if ((x==1.)&&(m==l)) return 1.;
  if ((x==1.)&&(m!=l)) return 0.;
  if ((x==-1.)&&(-m==l)) return powm1(n-l);
  if ((x==-1.)&&(-m!=l)) return 0.;
  if (abs(m)>n || abs(l)>n) return 0.;
  n0 = max(abs(m),abs(l));
  //if (theta < 0.) return coefs[l+N][m+N].t[n-n0];
  return coefs[m+N][l+N].t[n-n0];
}

WignerFunctionsC::WignerFunctionsC() {
  N = 0;
  coefs = NULL;
}

WignerFunctionsC::~WignerFunctionsC() {
  if (coefs != NULL) {
    for (int i=0 ; i<N+N+1 ; i++) delete [] (coefs[i]);
    delete [] (coefs);
  }
}

WignerFunctionsC::WignerFunctionsC(const int &nn, const Complex &X) {
  N = nn;
  x = X;
  int n,m,l,row,col;
  Complex d,dn,dnm1;
  double ee;
  coefs = new Vector*[N+N+1];
  for (int r=0 ; r<N+N+1 ; r++) coefs[r] = new Vector[N+N+1];
  row = col = 0;
  double lns[N+N+1];
  lns[0] = lns[1] = 0.;
  for (int i=1 ; i<N+N+1 ; i++) lns[i+1] = lns[i] + log(double(i+1));

  for (m=-N ; m<=N ; m++) {
    col = 0;
    for (l=-N ; l<=N ; l++) {
      n0 = max(abs(m),abs(l));
      coefs[row][col] = Vector(N-n0+1);
      ee = l<m ? powm1(m-l) : 1.;
      coefs[row][col].t[0] = d = ee*pow(2.,-double(n0)) * sqrt(exp(lns[n0+n0]-lns[int(abs(m-l))]-lns[int(abs(m+l))])) * pow(1.-x,double(abs(m-l))*0.5)*pow(1.+x,double(abs(m+l))*0.5);
      dn = d;
      dnm1 = 0.;
      n = 1;
      for (n0 = n0 ; n0<N ; n0++) {
        if (n0==0) coefs[row][col].t[n] = d = double(n0+n0+1)*x*dn - dnm1;
        else d = 1./(double(n0) * sqrt(double((n0+1)*(n0+1))-double(m*m)) * sqrt(double((n0+1)*(n0+1))-double(l*l)))
          * (double(n0+n0+1)*(double(n0*n0+n0)*x-double(m*l))*dn - double(n0+1)*sqrt(double(n0*n0-m*m))*sqrt(double(n0*n0-l*l))*dnm1);
        dnm1 = dn;
        coefs[row][col].t[n] = dn = d;
        n++;
      }
      col++;
    }
    row++;
  }
}


void WignerFunctionsC::assign(const int &nn, const Complex &X) {
  if (coefs != NULL) {
    for (int i=0 ; i<N+N+1 ; i++) delete [] (coefs[i]);
    delete [] (coefs);
  }
  N = nn;
  x = X;
  int n,m,l,row,col;
  Complex d,dn,dnm1;
  double ee;
  coefs = new Vector*[N+N+1];
  for (int r=0 ; r<N+N+1 ; r++) coefs[r] = new Vector[N+N+1];
  row = col = 0;
  double lns[N+N+1];
  lns[0] = lns[1] = 0.;
  for (int i=1 ; i<N+N+1 ; i++) lns[i+1] = lns[i] + log(double(i+1));
  for (m=-N ; m<=N ; m++) {
    col = 0;
    for (l=-N ; l<=N ; l++) {
      n0 = max(abs(m),abs(l));
      coefs[row][col] = Vector(N-n0+1,0.);
      ee = l<m ? powm1(m-l) : 1.;
      coefs[row][col].t[0] = dn = ee*pow(2.,-double(n0)) * sqrt(exp(lns[n0+n0]-lns[int(abs(m-l))]-lns[int(abs(m+l))])) * pow(1.-x,double(abs(m-l))*0.5)*pow(1.+x,double(abs(m+l))*0.5);
      dnm1 = 0.;
      n = 1;
      for (n0 = n0 ; n0<N ; n0++) {
        if (n0==0) coefs[row][col].t[n] = d = double(n0+n0+1)*x*dn - dnm1;
        else coefs[row][col].t[n] = d = 1./(double(n0) * sqrt(double((n0+1)*(n0+1))-double(m*m)) * sqrt(double((n0+1)*(n0+1))-double(l*l)))
          * (double(n0+n0+1)*(double(n0*n0+n0)*x-double(m*l))*dn - double(n0+1)*sqrt(double(n0*n0-m*m))*sqrt(double(n0*n0-l*l))*dnm1);
        dnm1 = dn;
        dn = d;
        n++;
      }
      col++;
    }
    row++;
  }
}

Complex WignerFunctionsC::dlmn(const int &l, const int &m, const int &n) {
  if ((x==1.)&&(m==l)) return 1.;
  if ((x==1.)&&(m!=l)) return 0.;
  if ((x==-1.)&&(-m==l)) return powm1(n-l);
  if ((x==-1.)&&(-m!=l)) return 0.;
  if (abs(m)>n || abs(l)>n) return 0.;
  n0 = max(abs(m),abs(l));
  //if (theta < 0.) return coefs[l+N][m+N].t[n-n0];
  //cout << m+N << "\t" << l+N << "\t" << n-n0 << endl;
  return coefs[m+N][l+N].t[n-n0];
}
