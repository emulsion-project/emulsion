#include "Integration.h"
#include <iostream>
using namespace std;

RMatrix gaussLegQuad(const int &n,const double &x1, const double &x2) {
  RMatrix x(n,2);
  const double eps = 1e-14;
  double z1,z,xm,xl,pp,p3,p2,p1;
  int m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for(int i=0 ; i<=m ; i++) {
    z = cos(piDbl*(i+0.75)/(n+0.5));
    do {
      p1 = 1.;
      p2 = 0.;
      for(int j=0 ; j<n ; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.*j+1.)*z*p2-j*p3)/double(j+1);
      }
      pp = double(n)*(z*p1-p2)/(z*z-1.);
      z1 = z;
      z = z1-p1/pp;
    } while (fabs(z-z1)>eps);
    x.t[i][0]=xm-xl*z;
    x.t[n-i-1][0] = xm+xl*z;
    x.t[i][1] = 2.0*xl/((1.0-z*z)*pp*pp);
    x.t[n-i-1][1] = x.t[i][1];
  }
  return x;
}

Matrix gaussLegQuad(const int &n,const Complex &x1, const Complex &x2) {
  Matrix x(n,2);
  const double eps = 1e-14;
  Complex z1,z,xm,xl,pp,p3,p2,p1;
  int m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for(int i=0 ; i<=m ; i++) {
    z = cos(piDbl*(i+0.75)/(n+0.5));
    do {
      p1 = 1.;
      p2 = 0.;
      for(int j=0 ; j<n ; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.*j+1.)*z*p2-j*p3)/double(j+1);
      }
      pp = double(n)*(z*p1-p2)/(z*z-1.);
      z1 = z;
      z = z1-p1/pp;
      //cout << abs(z-z1) << endl;
    } while (abs(z-z1)>eps);
    x.t[i][0]=xm-xl*z;
    x.t[n-i-1][0] = xm+xl*z;
    x.t[i][1] = 2.0*xl/((1.0-z*z)*pp*pp);
    x.t[n-i-1][1] = x.t[i][1];
  }
  return x;
}

RMatrix gaussLagQuad(const int &n,const double &alf) {
  const int maxit = 40;
  const double eps = 1.0e-14;
  int ai,i,its,j;
  double p1,p2,p3,pp,z,z1;
  RMatrix x(n,2);
  for (i=0 ; i<n ; i++) {
    if (i == 0) z = (1.+alf)*(3.+0.92*alf)/(1.+2.4*double(n)+1.8*alf);
    else if (i == 1) z += (15.+6.25*alf)/(1.+0.9*alf+2.5*double(n));
    else {
      ai = i-1;
      z += ((1.+2.55*double(ai))/(1.9*double(ai))+1.26*double(ai)*alf/(1.+3.5*double(ai)))*(z-x.t[i-2][0])/(1.+0.3*alf);
    }
    for (its=0 ; its<maxit ; its++) {
      p1 = 1.;
      p2 = 0.;
      for (j=0 ; j<n ; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((double(j+j+1)+alf-z)*p2-(double(j)+alf)*p3)/double(j+1);
      }
      pp = (double(n)*p1-(double(n)+alf)*p2)/z;
      z1 = z;
      z = z1-p1/pp;
      if (abs(z-z1) <= eps) break;
    }
    x.t[i][0] = z;
    x.t[i][1] = -exp(gammaln(alf+double(n))-gammaln(double(n)))/(pp*double(n)*p2);
  }
  return x;
}
