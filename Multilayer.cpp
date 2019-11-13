#include "Multilayer.h"
#include <iostream>
using namespace std;

Multilayer::Multilayer() {}

Multilayer::Multilayer(const int &nbrMed, const RVector &positions, const IVector &meds, const double &lambd) {
  setParameters(nbrMed,positions,meds,lambd);
}

void Multilayer::setParameters(const int &nbrMed, const RVector &positions, const IVector &meds, const double &lambd) {
  lambda = lambd;
  nbrMediums = nbrMed;
  Zi = positions;
  mediums = meds;
  ni = ki = Vector(nbrMediums);
  ri = new RefractiveIndex[nbrMediums];
  for (int i=0 ; i<nbrMediums ; i++) {
    ri[i].setMaterial(mediums.t[i],1.);
    ni.t[i] = ri[i].calcIndex(lambda);
    ki.t[i] = 2.*piDbl*ni.t[i]/lambda;
  }
}

void Multilayer::setLambda(const double &lambd) {
  lambda = lambd;
  for (int i=0 ; i<nbrMediums ; i++) {
    ni.t[i] = ri[i].calcIndex(lambda);
    ki.t[i] = 2.*piDbl*ni.t[i]/lambda;
  }
  return;
}

void Multilayer::setKxy(const Complex &k) {
  kxy = k;
  return;
}

int Multilayer::getMedium(const double &z) {
  if (nbrMediums == 1) return 0;
  for (int i=0 ; i<nbrMediums-1 ; i++) if (z<=Zi.t[i]) return i;
  return nbrMediums-1;
}

Complex Multilayer::calcKz(const int &lay, const bool &dirP) {
  Complex K = sqrt(ki.t[lay]*ki.t[lay]-kxy*kxy);
  if (!dirP) K*=-1.;
  if ((dirP && K.im<0.) || (!dirP && K.im > 0.)) K *= -1.;
  return K;
}

Vector Multilayer::sourceExtOut(const Vector &a, const int &pol) {
  Vector aa = calcSMatrix(pol)*a;
  return aa;
}

Vector Multilayer::sourceExtIn(const Vector &a, const int &pol, const double &z) {
  Vector aPM(2,0.),aa;
  int lay = getMedium(z);
  if (lay == 0) {
    aa = calcSMatrix(pol)*a;
    kz = calcKz(0,true);
    aPM.t[0] = a.t[0]*exp(i_*kz*(z-Zi.t[0]));
    aPM.t[1] = aa.t[0]*exp(-i_*kz*(z-Zi.t[0]));
    return aPM;
  }
  else if (lay == nbrMediums-1) {
    aa = calcSMatrix(pol)*a;
    kz = calcKz(nbrMediums-1,true);
    aPM.t[0] = aa.t[1]*exp(i_*kz*(z-Zi.t[nbrMediums-2]));
    aPM.t[1] = a.t[1]*exp(-i_*kz*(z-Zi.t[nbrMediums-2]));
    return aPM;
  }
  Matrix Sl,Su;
  Sl = calcSMatrixL(pol,z);
  Su = calcSMatrixU(pol,z);
  Complex cst = 1.-Su.t[0][0]*Sl.t[1][1];
  aPM.t[0] = (Sl.t[1][0]*a.t[0]+Sl.t[1][1]*Su.t[0][1]*a.t[1])/cst;
  aPM.t[1] = (Su.t[0][0]*Sl.t[1][0]*a.t[0]+Su.t[0][1]*a.t[1])/cst;

  return aPM;
}

Vector Multilayer::sourceInt(const int &pol, const Vector &a, const double &zs, const double &z) { // a = aP & aM
  aS=Vector(2,0.);
  Vector aSi(2,0.);
  Vector aa(2,0.);
  int lay = getMedium(zs);
  Matrix Sl,Su;
  Matrix S;

  if (lay == 0) {
    aS.t[0] = 0.;
    kz = calcKz(0,true);
    aS.t[1] = calcSMatrix(pol).t[0][0]*a.t[0]*exp(2.*i_*kz*(Zi.t[0]-zs));
  }
  else if (lay == nbrMediums-1) {
    aS.t[1] = 0.;
    kz = calcKz(nbrMediums-1,true);
    aS.t[0] = calcSMatrix(pol).t[1][1]*a.t[1]*exp(2.*i_*kz*(zs-Zi.t[nbrMediums-2]));
  }
  else {
    Sl = calcSMatrixL(pol,zs);
    Su = calcSMatrixU(pol,zs);
    Complex cst = 1.-Su.t[0][0]*Sl.t[1][1];
    aS.t[0] = Sl.t[1][1]*(Su.t[0][0]*a.t[0]+a.t[1])/cst;
    aS.t[1] = Su.t[0][0]*(Sl.t[1][1]*a.t[1]+a.t[0])/cst;
  }
  int layi = getMedium(z);
  if (layi == lay) {
    kz = calcKz(lay,true);
    aSi.t[0] = aS.t[0]*exp(i_*kz*(z-zs));
    aSi.t[1] = aS.t[1]*exp(-i_*kz*(z-zs));
  }
  /*else if (layi == 0) {
    kz = calcKz(0,true);
    if (zs>Zi.t[nbrMediums-2]) {
      Complex kz0 = calcKz(nbrMediums-1,true);
      S = calcSMatrix(pol);
      aSi.t[1] = S.t[0][1]*(aS.t[1]+a.t[1])*exp(i_*kz0*(zs-Zi.t[nbrMediums-2]))*exp(-i_*kz*(z-Zi.t[0]));
    }
    else {
      S = calcSMatrixL(pol,zs);
      aSi.t[1] = S.t[0][1]*(aS.t[1]+a.t[1])*exp(-i_*kz*(z-Zi.t[0]));
    }
    aSi.t[0] = 0.;
  }*/
  /*else if (z > Zi.t[nbrMediums-2]) {
    kz = calcKz(nbrMediums-1,true);
    S = calcSMatrixU(pol,zs);
    aSi.t[0] = S.t[1][0]*(aS.t[0]+a.t[0])*exp(i_*kz*(z-Zi.t[nbrMediums-2]));
    aSi.t[1] = 0.;
  }*/
  else if (lay == 0) {
    kz = calcKz(0,true);
    aa.t[0] = a.t[0]*exp(i_*kz*(Zi.t[0]-zs));
    aa.t[1] = 0.;
    aSi = sourceExtIn(aa,pol,z);
  }
  else if (lay == nbrMediums-1) {
    kz = calcKz(nbrMediums-1,true);
    aa.t[1] = a.t[1]*exp(i_*kz*(zs-Zi.t[nbrMediums-2]));
    aa.t[0] = 0.;
    aSi = sourceExtIn(aa,pol,z);
  }
  else if (z<zs) {
    S = calcSMatrixInt(pol,z,zs);
    if (S.t[1][0] < 1e-299) S.t[1][0] = 1e-299;
    aSi.t[0] = (aS.t[0]-S.t[1][1]*(aS.t[1]+a.t[1]))/S.t[1][0];
    aSi.t[1] = S.t[0][0]*aSi.t[0] + S.t[0][1]*(aS.t[1]+a.t[1]);
  }
  else if (z>zs) {
    S = calcSMatrixInt(pol,zs,z);
    if (S.t[0][1] < 1e-299) S.t[0][1] = 1e-299;
    aSi.t[1] = (aS.t[1]-S.t[0][0]*(aS.t[0]+a.t[0]))/S.t[0][1];
    aSi.t[0] = S.t[1][0]*(aS.t[0]+a.t[0])+S.t[1][1]*aSi.t[1];
  }
  return aSi;
}

Vector Multilayer::calcFresnelCoefsTE(const int &i, const int &j) {
  Vector coefs(2,0.);
  if (i<j) {
    Kzi = calcKz(i,true);
    Kzj = calcKz(j,true);
  }
  else {
    Kzi = calcKz(i,false);
    Kzj = calcKz(j,false);
  }
  coefs.t[0] = (Kzi - Kzj)/(Kzj + Kzi); //rTE, ra, rs
  coefs.t[1] = (2.*Kzi)/(Kzj + Kzi); // tTE
  return coefs;
}

Vector Multilayer::calcFresnelCoefsTM(const int &i, const int &j) {
  Complex Kzi,Kzj;
  Vector coefs(2,0.);
  if (i<j) {
    Kzi = calcKz(i,true);
    Kzj = calcKz(j,true);
  }
  else {
    Kzi = calcKz(i,false);
    Kzj = calcKz(j,false);
  }
  coefs.t[0] = -(ni.t[i]*ni.t[i]*Kzj - ni.t[j]*ni.t[j]*Kzi)/(ni.t[j]*ni.t[j]*Kzi + ni.t[i]*ni.t[i]*Kzj); //rTM, rb, rp
  coefs.t[1] = (2.*ni.t[i]*ni.t[j]*Kzi)/(ni.t[j]*ni.t[j]*Kzi + ni.t[i]*ni.t[i]*Kzj); // tTM
  return coefs;
}

Matrix Multilayer::calcSl(const int &l,const int &u,const int &pol) {
  Matrix S(2,2);
  Vector coefsUL, coefsLU;
  if (pol == 0) {
    coefsUL = calcFresnelCoefsTE(u,l);
    coefsLU = calcFresnelCoefsTE(l,u);
  }
  else {
    coefsUL = calcFresnelCoefsTM(u,l);
    coefsLU = calcFresnelCoefsTM(l,u);
  }
  S.t[0][0] = coefsLU.t[0];
  S.t[0][1] = coefsUL.t[1];
  S.t[1][0] = coefsLU.t[1];
  S.t[1][1] = coefsUL.t[0];
  return S;
}

Matrix Multilayer::calcSh(const int &i, const double &h) {
  Matrix S(2,2);
  kz = calcKz(i,true);
  S.t[0][0] = S.t[1][1] = 0.;
  S.t[0][1] = S.t[1][0] = exp(i_*kz*h);
  return S;
}

Matrix Multilayer::calcSMatrix(const int &pol) { // 0:TE 1:TM
  Matrix S,Su,Sl;
  Complex cst;
  S = Su = Sl = Matrix(2,2,0.);
  Su = calcSl(nbrMediums-2,nbrMediums-1,pol);
  if (nbrMediums == 1) return calcSh(0,0.);
  if (nbrMediums == 2) return calcSl(0,1,pol);

  for (int i=nbrMediums-2 ; i>0 ; i--) {
    Sl = calcSh(i,Zi.t[i]-Zi.t[i-1]);
    cst = 1./(1.-Su.t[0][0]*Sl.t[1][1]);
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]*Sl.t[1][0]*cst;
    S.t[0][1] = Sl.t[0][1]*Su.t[0][1]*cst;
    S.t[1][0] = Sl.t[1][0]*Su.t[1][0]*cst;
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]*Su.t[0][1]*cst;

    Su = S;
    Sl = calcSl(i-1,i,pol);

    cst = 1./(1.-Su.t[0][0]*Sl.t[1][1]);
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]*Sl.t[1][0]*cst;
    S.t[0][1] = Sl.t[0][1]*Su.t[0][1]*cst;
    S.t[1][0] = Sl.t[1][0]*Su.t[1][0]*cst;
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]*Su.t[0][1]*cst;

    Su = S;
  }
  return S;
}

Matrix Multilayer::calcSMatrixU(const int &pol,const double &z) { // 0:TE 1:TM
  Matrix S,Su,Sl;
  S = Su = Sl = Matrix(2,2,0.);
  int lay = getMedium(z);
  if (lay == 0) {
    Sl = calcSh(0,fabs(Zi.t[0]-z));
    Su = calcSMatrix(pol);
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    return S;
  }
  //if (lay == nbrMediums-1) return S;
  if (lay == nbrMediums-1) return calcSh(nbrMediums-1,fabs(z-Zi.t[nbrMediums-2]));

  Su = calcSl(nbrMediums-2,nbrMediums-1,pol);

  for (int i=nbrMediums-2 ; i>0 ; i--) {
    if (lay == i) Sl = calcSh(i,fabs(Zi.t[i]-z));
    else Sl = calcSh(i,fabs(Zi.t[i-1]-Zi.t[i]));

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    if (lay == i) break;

    Su = S;
    Sl = calcSl(i-1,i,pol);

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];

    Su = S;
  }
  return S;
}

Matrix Multilayer::calcSMatrixL(const int &pol,const double &z) { // 0:TE 1:TM
  Matrix S,Su,Sl;
  S = Su = Sl = Matrix(2,2,0.);
  int lay = getMedium(z);
  if (lay == nbrMediums-1) {
    Su = calcSh(nbrMediums-1,fabs(z-Zi.t[nbrMediums-2]));
    Sl = calcSMatrix(pol);
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    return S;
  }
  //if (lay == 0) return S;
  if (lay == 0) return calcSh(0,fabs(z-Zi.t[0]));
  Sl = calcSl(0,1,pol);

  for (int i=1 ; i<nbrMediums ; i++) {
    if (lay == i) Su = calcSh(i,fabs(z-Zi.t[i-1]));
    else Su = calcSh(i,fabs(Zi.t[i]-Zi.t[i-1]));

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    if (lay == i) break;

    Sl = S;
    Su = calcSl(i,i+1,pol);

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];

    Sl = S;
  }
  return S;
}

Matrix Multilayer::calcSMatrixInt(const int &pol,const double &z1, const double &z2) { // 0:TE 1:TM
  Matrix S,Su,Sl;
  S = Su = Sl = Matrix(2,2,0.);
  int lay1 = getMedium(z1);
  int lay2 = getMedium(z2);
  if (lay1 == lay2) return calcSh(lay1,fabs(z2-z1));
  if (lay1 == 0) {
    Su = calcSMatrixL(pol,z2);
    Sl = calcSh(0,fabs(Zi.t[0]-z1));
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    return S;
  }

  if (lay2 == nbrMediums-1) {
    Sl = calcSMatrixU(pol,z1);
    Su = calcSh(nbrMediums-1,fabs(z2-Zi.t[nbrMediums-2]));
    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    return S;
  }
  Sl = calcSh(lay1,fabs(Zi.t[lay1]-z1));

  for (int i=lay1+1 ; i<nbrMediums ; i++) {
    Su = calcSl(i-1,i,pol);

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];

    Sl = S;

    if (lay2 == i) Su = calcSh(i,fabs(z2-Zi.t[i-1]));
    else Su = calcSh(i,fabs(Zi.t[i]-Zi.t[i-1]));

    S.t[0][0] = Sl.t[0][0] + Sl.t[0][1]*Su.t[0][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Sl.t[1][0];
    S.t[0][1] = Sl.t[0][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];
    S.t[1][0] = Sl.t[1][0]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[1][0];
    S.t[1][1] = Su.t[1][1] + Su.t[1][0]*Sl.t[1][1]/(1.-Su.t[0][0]*Sl.t[1][1])*Su.t[0][1];

    if (lay2 == i) break;
    Sl = S;
  }
  return S;
}
