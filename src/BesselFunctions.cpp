#include "BesselFunctions.h"
#include <iostream>
using namespace std;

BesselFunctions::BesselFunctions() {}

/** Function computing the Bessel functions of first kind with real positive order nu = N+alpha and complex argument zz **/
void BesselFunctions::calcCJ(const int &N, const double &alpha, const Complex &zz) {
  z = zz;
  if (z.im<0.) z.im *= -1.;
  /// Initialization of the vector containing the values of Jnu
  CJ = Vector(N+1,0.);
  /// Calculation of the backward recurrence starting order eta
  Complex j1=0.,j2=1.,j3=0.;
  n = N;
  while (abs(j3)<2e15) {
    j3 = double(n+n-2)/z*j2 - j1;
    j1 = j2; j2 = j3; n++;
  }
  eta = n;

  Vector s(eta+1),r(eta+1),l(eta+1),lambda(eta+1);
  lambda.t[1] = -2.*i_*(alpha+1.);
  l.t[1] = 1.;
  for (n=1 ; n<=eta-1 ; n++) {
    l.t[n+1] = (double(n)+2.*alpha)/double(n+1)*l.t[n];
    lambda.t[n+1] = 2.*powi(-n-1)*(alpha+double(n+1))*l.t[n+1];
  }

  s.t[eta] = r.t[eta] = 0.;
  for (n=eta ; n>=1 ; n--) {
    r.t[n-1] = 1./(2.*(alpha+double(n))/z - r.t[n]);
    s.t[n-1] = r.t[n-1]*(lambda.t[n] + s.t[n]);
  }
  Complex ss = pow(z*0.5,alpha)*exp(-i_*z)/exp(gammaln(1.+alpha));
  CJ.t[0] = ss/(1.+s.t[0]);
  for (n=1 ; n<=N ; n++) CJ.t[n] = r.t[n-1]*CJ.t[n-1];

  if (zz.im < 0.) for (n=0 ; n<=N ; n++) CJ.t[n].im *= -1.;
  return;
}

/** Function computing the Bessel functions of first kind with integer positive order N and complex argument zz **/
void BesselFunctions::calcCJ(const int &N, const Complex &zz) {
  z = zz;
  if (z.im<0.) z.im *= -1.;
  /// Initialization of the vector containing the values of Jnu
  CJ = Vector(N+1,0.);
  /// Calculation of the backward recurrence starting order eta
  Complex j1=0.,j2=1.,j3=0.;
  n = N;
  while (abs(j3)<2e15) {
    j3 = double(n+n-2)/z*j2 - j1;
    j1 = j2; j2 = j3; n++;
  }
  eta = n;

  Vector s(eta+1),r(eta+1),lambda(eta+1);
  for (n=1 ; n<=eta-1 ; n++) lambda.t[n] = 2.*powi(-n);

  s.t[eta] = r.t[eta] = 0.;
  for (n=eta ; n>=1 ; n--) {
    r.t[n-1] = 1./(2.*double(n)/z - r.t[n]);
    s.t[n-1] = r.t[n-1]*(lambda.t[n] + s.t[n]);
  }
  CJ.t[0] = exp(-i_*z)/(1.+s.t[0]);
  for (n=1 ; n<=N ; n++) CJ.t[n] = r.t[n-1]*CJ.t[n-1];

  if (zz.im < 0.) for (n=0 ; n<=N ; n++) CJ.t[n].im *= -1.;
  return;
}

void BesselFunctions::calcJ(const int &N, const double &alpha, const double &xx) {
  x = xx;
  if (x<0.) cout << "invalid argument" << endl;
  /// Initialization of the vector containing the values of Jnu
  J = RVector(N+1,0.);
  /// Calculation of the backward recurrence starting order eta
  double j1=0.,j2=1.,j3=0.;
  n = N;
  while (abs(j3)<2e15) {
    j3 = double(n+n-2)/x*j2 - j1;
    j1 = j2; j2 = j3; n++;
  }
  eta = n;
  if (eta%2 != 0) eta++;

  RVector s(eta+1),r(eta+1),l(eta+1),lambda(eta+1);
  lambda.t[2] = alpha+2.;
  lambda.t[1] = 0.;
  l.t[1] = 1.;
  for (n=1 ; n<eta/2 ; n++) {
    l.t[n+1] = (double(n)+alpha)/double(n+1)*l.t[n];
    lambda.t[2*n+2] = (alpha+2.*double(n+1))*l.t[n+1];
    lambda.t[2*n+1] = 0.;
  }

  s.t[eta] = r.t[eta] = 0.;
  for (n=eta ; n>=1 ; n--) {
    r.t[n-1] = 1./(2.*(alpha+double(n))/x - r.t[n]);
    s.t[n-1] = r.t[n-1]*(lambda.t[n] + s.t[n]);
  }
  double ss = pow(x*0.5,alpha)/exp(gammaln(1.+alpha));
  J.t[0] = ss/(1.+s.t[0]);
  for (n=1 ; n<=N ; n++) J.t[n] = r.t[n-1]*J.t[n-1];

  //if (xx < 0.) for (n=0 ; n<=N ; n++) J.t[n] *= exp(i_*piDbl*x);
  return;
}

void BesselFunctions::calcJ(const int &N, const double &xx) {
  x = xx;
  if (x<0.) x *= -1.;
  /// Initialization of the vector containing the values of Jnu
  J = RVector(N+1,0.);
  /// Calculation of the backward recurrence starting order eta
  double j1=0.,j2=1.,j3=0.;
  n = N;
  while (abs(j3)<2e15) {
    j3 = double(n+n-2)/x*j2 - j1;
    j1 = j2; j2 = j3; n++;
  }
  eta = n;
  if (eta%2 != 0) eta++;

  RVector s(eta+1),r(eta+1),lambda(eta+1,0.);
  for (n=0 ; n<eta/2 ; n++) lambda.t[2*n+2] = 2.;

  s.t[eta] = r.t[eta] = 0.;
  for (n=eta ; n>=1 ; n--) {
    r.t[n-1] = 1./(2.*double(n)/x - r.t[n]);
    s.t[n-1] = r.t[n-1]*(lambda.t[n] + s.t[n]);
  }
  J.t[0] = 1./(1.+s.t[0]);
  for (n=1 ; n<=N ; n++) J.t[n] = r.t[n-1]*J.t[n-1];

  if (xx < 0.) for (n=0 ; n<=N ; n++) J.t[n] *= powm1(n);
  return;
}

SphericalBessel::SphericalBessel() {}

SphericalBessel::SphericalBessel(const int &NN, const Complex &zz) {
  J = Vector(NN+1);
  dJ = Vector(NN+1);
  Vector R(NN+1);
  Complex r;
  int NP = 300;
  double eps = 1e-13;
  Complex z;
  if (abs(zz)<eps) {
    J.t[0] = 1.;
    for (int i=1 ; i<=NN ; i++) J.t[i] = dJ.t[i] = 0.;
  }
  else {
    z = zz;
    J.t[0] = sin(z)/z;
    r = 0.;
    for (int i=NP ; i>0 ; i--) {
      r = 1./((2.*double(i)-1.)/z -r);
      if (i<=NN+1) R.t[i-1] = r;
    }
    dJ.t[0] = cos(z);
    for (int i=1 ; i<=NN ; i++) {
      J.t[i] = J.t[i-1]*R.t[i];
      dJ.t[i] = z*J.t[i-1] - double(i)*J.t[i];
    }
  }
  Y = Vector(NN+1);
  dY = Vector(NN+1);
  Y.t[0] = -cos(z)/z;
  dY.t[0] = sin(z);
  if (NN>0) {
    Y.t[1] = Y.t[0]/z - sin(z)/z;
    dY.t[1] = z*Y.t[0] - Y.t[1];
  }
  for (int i=2 ; i<=NN ; i++) {
    Y.t[i] = double(i+i-1)*Y.t[i-1]/z - Y.t[i-2];
    dY.t[i] = z*Y.t[i-1] - double(i)*Y.t[i];
  }
}

void SphericalBessel::assign(const int &NN, const Complex &zz) {
  J = Vector(NN+1);
  dJ = Vector(NN+1);
  Vector R(NN+1);
  Complex r;
  int NP = 500;
  double eps = 1e-13;
  Complex z;
  if (abs(zz)<eps) {
    J.t[0] = 1.;
    for (int i=1 ; i<=NN ; i++)
      J.t[i] = dJ.t[i] = 0.;
  }
  else {
    z = zz;
    J.t[0] = sin(z)/z;
    r = 0.;
    for (int i=NP ; i>0 ; i--) {
      r = 1./((2.*double(i)-1.)/z -r);
      if (i<=NN+1) R.t[i-1] = r;
    }
    dJ.t[0] = cos(z);
    for (int i=1 ; i<=NN ; i++) {
      J.t[i] = J.t[i-1]*R.t[i];
      dJ.t[i] = z*J.t[i-1] - double(i)*J.t[i];
    }
  }
  Y = Vector(NN+1);
  dY = Vector(NN+1);
  Y.t[0] = -cos(z)/z;
  dY.t[0] = sin(z);
  if (NN>0) {
    Y.t[1] = Y.t[0]/z - sin(z)/z;
    dY.t[1] = z*Y.t[0] - Y.t[1];
  }
  for (int i=2 ; i<=NN ; i++) {
    Y.t[i] = double(i+i-1)*Y.t[i-1]/z - Y.t[i-2];
    dY.t[i] = z*Y.t[i-1] - double(i)*Y.t[i];
  }
}

Complex SphericalBessel::jn(const int &n) { return J.t[n]; }
Complex SphericalBessel::yn(const int &n) { return Y.t[n]; }
Complex SphericalBessel::h1n(const int &n) { return J.t[n] + i_*Y.t[n]; }
Complex SphericalBessel::h2n(const int &n) { return J.t[n] - i_*Y.t[n]; }
Complex SphericalBessel::djn(const int &n) { return dJ.t[n]; }
Complex SphericalBessel::dyn(const int &n) { return dY.t[n]; }
Complex SphericalBessel::dh1n(const int &n) { return dJ.t[n] + i_*dY.t[n]; }
Complex SphericalBessel::dh2n(const int &n) { return dJ.t[n] - i_*dY.t[n]; }

Bessel::Bessel() {}

Bessel::Bessel(const int &NN, const Complex &z) {
  J = Vector(NN+1,0.);
  Y = Vector(NN+1,0.);
  Complex ll=0.,j1,j2,j3;
  double lim = 1e15;
  j1 = 0.;
  j2 = 1.;
  int n = NN;
  while (abs(ll)<lim) {
    j3 = 2.*double(n-1)/z*j2 - j1;
    ll = j3;
    j1 = j2;
    j2 = j3;
    n++;
  }
  int nn = n+1;
  Vector j(nn+1,0.);
  j.t[nn] = 0.;
  j.t[nn-1] = 1.;
  ll=0.;
  for (n=nn-2 ; n>=0 ; n--) {
    j.t[n] = 2.*double(n+1)/z*j.t[n+1] - j.t[n+2];
    if (n%2 == 0) ll += j.t[n];
  }
  for (n=0 ; n<=NN ; n++) J.t[n] = j.t[n]/(2.*ll-j.t[0]);
  for (n=nn-1 ; n>=0 ; n--) j.t[n] /= 2.*ll-j.t[0];

  Y.t[0] = 2./piDbl*(log(0.5*z)+euler)*J.t[0];
  for (n=1 ; n<=nn/2-1 ; n++) Y.t[0] -= 4./piDbl*powm1(n)*j.t[n*2]/double(n);
  for (n=1 ; n<=NN ; n++) Y.t[n] = (Y.t[n-1]*J.t[n]-2./(piDbl*z))/J.t[n-1];
}

void Bessel::assign(const int &NN, const Complex &z) {
  J = Vector(NN+1,0.);
  if (abs(z) == 0.) {
    J.t[0] = 1.;
    return;
  }
  //Y = Vector(NN+1,0.);
  Complex ll=0.,j1,j2,j3;
  double lim = 2e15;
  j1 = 0.;
  j2 = 1.;
  int n = 1;
  while (abs(ll)<lim) {
    j3 = 2.*double(n-1)/z*j2 - j1;
    ll = j3;
    j1 = j2;
    j2 = j3;
    n++;
  }
  int nn = n-1;
  Vector j(nn+1,0.);
  j.t[nn] = 0.;
  j.t[nn-1] = 1.;
  ll=0.;
  for (n=nn-2 ; n>=0 ; n--) {
    j.t[n] = 2.*double(n+1)/z*j.t[n+1] - j.t[n+2];
    if (n%2 == 0) ll += j.t[n];
  }
  for (n=0 ; n<=NN ; n++) J.t[n] = j.t[n]/(2.*ll-j.t[0]);
  for (n=nn-1 ; n>=0 ; n--) j.t[n] /= 2.*ll-j.t[0];

  /*Y.t[0] = 2./piDbl*(log(0.5*z)+euler)*J.t[0];
  for (n=1 ; n<=nn/2-1 ; n++) Y.t[0] -= 4./piDbl*powm1(n)*j.t[n*2]/double(n);
  for (n=1 ; n<=NN ; n++) Y.t[n] = (Y.t[n-1]*J.t[n]-2./(piDbl*z))/J.t[n-1];*/
  return;
}

Complex Bessel::Jn(const int &n) { return J.t[n]; }

Complex besselJn(const int &n, const Complex &z) {
  RMatrix wx;
  Complex jn,jnn=piDbl;
  int i=200;
  wx = gaussLegQuad(i,0.,piDbl);
  jn=0.;
  for (int j=0 ; j<i ; j++) jn+= wx.t[j][1]*cos(z*sin(wx.t[j][0])-double(n)*wx.t[j][0]);
  jn /= piDbl;
  jnn = jn;
  return jn;
}

double besselI0(const double &x) {
  float ax,ans;
  double y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
  }
  return ans;
}

double besselK0(const double & x) {
  double y,ans;
  if (x <= 2.0) {
    y=x*x/4.0;
    ans=(-log(x/2.0)*besselI0(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5))))));
  }
  else {
    y=2.0/x;
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))));
  }
  return ans;
}

double besselI1(const double &x) {
  double ax,ans;
  double y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  else {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}


double besselK1(const double &x) {
  double y,ans;
  if (x <= 2.0) {
    y=x*x/4.0;
    ans=(log(x/2.0)*besselI1(x))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))));
  }
  else {
    y=2.0/x;
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619 +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
  return ans;
}


double besselKn(const int &nn, const double &x) {
  int n = abs(nn);
  int j;
  double bk,bkm,bkp,tox;
  if (n == 0) return besselK0(x);
  if (n == 1) return besselK1(x);
  tox=2./x;
  bkm=besselK0(x);
  bk=besselK1(x);
  for (j=1;j<n;j++) {
    bkp=bkm+j*tox*bk;
    bkm=bk;
    bk=bkp;
  }
  return bk;
}


double besselIn(int n, float x) {
  double ACC = 40.0; // Make larger to increase accuracy.
  double BIGNO = 1.0e10;
  double BIGNI = 1.0e-10;
  int j;
  double bi,bim,bip,tox,ans;
  if (n == 0) return besselI0(x);
  if (n == 1) return besselI1(x);
  if (x == 0.0) return 0.0;
  else {
    tox=2.0/fabs(x);
    bip=ans=0.0;
    bi=1.0;
    for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
      bim=bip+j*tox*bi;
      bip=bi;
      bi=bim;
      if (fabs(bi) > BIGNO) {
        ans *= BIGNI;
        bi *= BIGNI;
        bip *= BIGNI;
      }
      if (j == n) ans=bip;
    }
    ans *= besselI0(x)/bi;
    return x < 0.0 && (n & 1) ? -ans : ans;
  }
}
