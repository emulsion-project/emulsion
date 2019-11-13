#include "Multilayer.h"
#include <iostream>
using namespace std;

Multilayer::Multilayer() {}

Multilayer::Multilayer(const string &fname, const double &lambd) {
  setParameters(fname,lambd);
}

void Multilayer::setParameters(const string &fname, const double &lambd) {
  string str,fileName;
  ifstream ifs(fname.c_str());
  lambda = lambd;
  ifs >> str >> nbrMediums >> str;
  int i;

  if (nbrMediums>1) {
    Zi = RVector(nbrMediums-1);
    for (i=0 ; i<nbrMediums-1 ; i++) ifs >> Zi.t[i];
  }
  else {
    Zi = RVector(1);
    ifs >> Zi.t[0];
  }
  ifs >> str;
  ni = ki = Vector(nbrMediums);
  ri = new RefractiveIndex[nbrMediums];
  mediums = new int[nbrMediums];
  for (i=0 ; i<nbrMediums ; i++) {
    ifs >> mediums[i];
    ri[i].setMaterial(mediums[i],1.);
    ni.t[i] = ri[i].calcIndex(lambda);
    ki.t[i] = 2.*piDbl*ni.t[i]/lambda;
  }
}

int Multilayer::getMedium(const double &z) {
  if (nbrMediums == 1) return 0;
  for (int i=0 ; i<nbrMediums-1 ; i++) if (z<=Zi.t[i]) return i;
  return nbrMediums-1;
}

Complex Multilayer::calcKz(const int &lay, const bool &dirP) {
  Complex K = sqrt(ki.t[lay]*ki.t[lay]-kx*kx-ky*ky);
  if (!dirP) K*=-1.;
  if ((dirP && K.im<0.) || (!dirP && K.im > 0.)) K *= -1.;
  return K;
}

Complex Multilayer::getAngle(const int &l1, const Complex &beta1, const double &alpha1, const int &l2) {
  kx = ki.t[l1]*cos(alpha1)*sin(beta1);
  ky = ki.t[l1]*sin(alpha1)*sin(beta1);
  if (beta1 <= piDblby2) return acos(calcKz(l2,true)/ki.t[l2]);
  else return acos(calcKz(l2,false)/ki.t[l2]);
}

void Multilayer::setSourceExt(const Complex &e0, const Complex &beta, const double &alpha) {
  a0 = Vector(2,0.);
  alpha0 = alpha;
  beta0 = beta;
  if (beta0<piDbl/2.) {
    kx = ki.t[0]*cos(alpha0)*sin(beta0);
    ky = ki.t[0]*sin(alpha0)*sin(beta0);
    a0.t[0] = e0; // a1+
  }
  else {
    kx = ki.t[nbrMediums-1]*cos(alpha0)*sin(beta0);
    ky = ki.t[nbrMediums-1]*sin(alpha0)*sin(beta0);
    a0.t[1] = e0; //a2-
  }
  return;
}

Vector Multilayer::calcIncAmplOut(const int &pol) {
  return calcSMatrix(pol)*a0;
}

Vector3d Multilayer::sourceExtFields(const Complex &e0, const Complex &beta, const double &alpha,const int &pol, const Vector3d &pos) {
  a0 = Vector(2,0.);
  alpha0 = alpha;
  beta0 = beta;
  if (beta0<piDblby2) {
    kx = ki.t[0]*cos(alpha0)*sin(beta0);
    ky = ki.t[0]*sin(alpha0)*sin(beta0);
    a0.t[0] = e0; // a1+
  }
  else {
    kx = ki.t[nbrMediums-1]*cos(alpha0)*sin(beta0);
    ky = ki.t[nbrMediums-1]*sin(alpha0)*sin(beta0);
    a0.t[1] = e0; //a2-
  }

  double z = real(pos.t[2]);
  Vector a = exp(i_*(kx*pos.t[0]+ky*pos.t[1]))*calcIncAmplIn(z,pol);
  int lay = getMedium(z);
  //cout << a0 << endl;
  kz = calcKz(lay,true);
  Complex bP = acos(kz/ki.t[lay]), bM = acos(-kz/ki.t[lay]);
  Vector3d EP(0.),EM(0.);

  EP.t[2-pol] = a.t[0];
  EM.t[2-pol] = a.t[1];
  EP = toCartesianCoordinates(EP,bP,alpha0);
  EM = toCartesianCoordinates(EM,bM,alpha0);

  return EP+EM;
}

Vector Multilayer::sourceExtAmpl(const Complex &e0, const Complex &beta, const double &alpha,const int &pol, const Vector3d &pos) {
  a0 = Vector(2,0.);
  alpha0 = alpha;
  beta0 = beta;
  if (beta0<piDblby2) {
    kx = ki.t[0]*cos(alpha0)*sin(beta0);
    ky = ki.t[0]*sin(alpha0)*sin(beta0);
    a0.t[0] = e0; // a1+
  }
  else {
    kx = ki.t[nbrMediums-1]*cos(alpha0)*sin(beta0);
    ky = ki.t[nbrMediums-1]*sin(alpha0)*sin(beta0);
    a0.t[1] = e0; //a2-
  }
  double z = real(pos.t[2]);
  Vector a = exp(i_*(kx*pos.t[0]+ky*pos.t[1]))*calcIncAmplIn(z,pol);
  return a;
}

Vector Multilayer::calcIncAmplIn(const double &z, const int &pol) {
  Vector aPM(2,0.),aa;
  int lay = getMedium(z);
  if (lay == 0) {
    aa = calcSMatrix(pol)*a0;
    kz = calcKz(0,true);
    aPM.t[0] = a0.t[0]*exp(i_*kz*(z-Zi.t[0]));
    aPM.t[1] = aa.t[0]*exp(-i_*kz*(z-Zi.t[0]));
    return aPM;
  }
  else if (lay == nbrMediums-1) {
    aa = calcSMatrix(pol)*a0;
    kz = calcKz(nbrMediums-1,true);
    aPM.t[0] = aa.t[1]*exp(i_*kz*(z-Zi.t[nbrMediums-2]));
    aPM.t[1] = a0.t[1]*exp(-i_*kz*(z-Zi.t[nbrMediums-2]));
    return aPM;
  }
  Matrix Sl,Su;
  aa = a0;
  Sl = calcSMatrixL(pol,z);
  Su = calcSMatrixU(pol,z);
  Complex cst = 1.-Su.t[0][0]*Sl.t[1][1];
  aPM.t[0] = (Sl.t[1][0]*aa.t[0]+Sl.t[1][1]*Su.t[0][1]*aa.t[1])/cst;
  aPM.t[1] = (Su.t[0][0]*Sl.t[1][0]*aa.t[0]+Su.t[0][1]*aa.t[1])/cst;

  return aPM;
}

Vector Multilayer::setSourceInt(const int &pol, const Complex &aP, const Complex &aM, const Complex &betaP, const double &alphaP, const double &z) {
  a=aS=Vector(2,0.);
  int lay = getMedium(z);
  kx = ki.t[lay]*cos(alphaP)*sin(betaP);
  ky = ki.t[lay]*sin(alphaP)*sin(betaP);
  a.t[0] = aP;
  a.t[1] = aM;

  if (lay == 0) {
    aS.t[0] = 0.;
    kz = calcKz(0,true);
    aS.t[1] = calcSMatrix(pol).t[0][0]*a.t[0]*exp(2.*i_*kz*(Zi.t[0]-z));
    //aS.t[1] = calcSMatrix(pol).t[0][0]*a.t[0];//*exp(2.*i_*kz*(Zi.t[0]-z));
  }
  else if (lay == nbrMediums-1) {
    aS.t[1] = 0.;
    kz = calcKz(nbrMediums-1,true);
    aS.t[0] = calcSMatrix(pol).t[1][1]*a.t[1]*exp(2.*i_*kz*(z-Zi.t[nbrMediums-2]));
  }
  else aS = calcSourceAmplIn(z,pol);
  return aS;
}

Vector Multilayer::setSourceIntAmpl(const int &pol, const Complex &aP, const Complex &aM, const Complex &betaP, const double &alphaP, const double &zs, const Vector3d &pos) {
  a=aS=Vector(2,0.);
  Vector aSi(2,0.);
  int lay = getMedium(zs);
  kx = ki.t[lay]*cos(alphaP)*sin(betaP);
  ky = ki.t[lay]*sin(alphaP)*sin(betaP);

  a.t[0] = aP;
  a.t[1] = aM;

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
  else aS = calcSourceAmplIn(zs,pol);

  Matrix S;
  double z = real(pos.t[2]);
  if (z<zs) S = calcSMatrixInt(pol,z,zs);
  else S = calcSMatrixInt(pol,zs,z);
  int layi = getMedium(z);
  if (layi == lay) {
    kz = calcKz(lay,true);
    aSi.t[0] = aS.t[0]*exp(i_*kz*(z-zs));
    aSi.t[1] = aS.t[1]*exp(-i_*kz*(z-zs));
  }
  else if (layi == 0) {
    kz = calcKz(0,true);
    aSi.t[1] = S.t[0][1]*(aS.t[1]+a.t[1])*exp(-i_*kz*(z-Zi.t[0]));
    aSi.t[0] = 0.;
  }
  else if (layi == nbrMediums-1) {
    S = calcSMatrixU(pol,zs);
    kz = calcKz(nbrMediums-1,true);
    aSi.t[0] = S.t[1][0]*(aS.t[0]+a.t[0])*exp(i_*kz*(z-Zi.t[nbrMediums-2]));
    aSi.t[1] = 0.;
  }
  else if (z<zs) {
    aSi.t[0] = (aS.t[0]-S.t[1][1]*(aS.t[1]+a.t[1]))/S.t[1][0];
    aSi.t[1] = (S.t[0][0]*aS.t[0]+(S.t[0][1]*S.t[1][0]-S.t[0][0]*S.t[1][1])*(aS.t[1]+a.t[1]))/S.t[1][0];
  }
  else if (z>zs) {
    aSi.t[0] = ((S.t[0][1]*S.t[1][0]-S.t[0][0]*S.t[1][1])*(aS.t[0]+a.t[0])+S.t[1][1]*aS.t[1])/S.t[0][1];
    aSi.t[1] = (aS.t[1]-S.t[0][0]*(aS.t[0]+a.t[0]))/S.t[0][1];
  }
  aSi = exp(i_*(kx*pos.t[0]+ky*pos.t[1]))*aSi;
  return aSi;
}

Vector3d Multilayer::setSourceIntFields(const int &pol, const Complex &aP, const Complex &aM, const Complex &betaP, const double &alphaP, const double &zs, const Vector3d &pos) {
  a=aS=Vector(2,0.);
  Vector aSi(2,0.);
  int lay = getMedium(zs);


  kx = ki.t[lay]*cos(alphaP)*sin(betaP);
  ky = ki.t[lay]*sin(alphaP)*sin(betaP);
  Complex kz;

  a.t[0] = aP;
  a.t[1] = aM;

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
  else aS = calcSourceAmplIn(zs,pol);

  Matrix S;
  double z = real(pos.t[2]);
  if (z<zs) S = calcSMatrixInt(pol,z,zs);
  else S = calcSMatrixInt(pol,zs,z);
  int layi = getMedium(z);
  if (layi == lay) {
    kz = calcKz(lay,true);
    aSi.t[0] = aS.t[0]*exp(i_*kz*(z-zs));
    aSi.t[1] = aS.t[1]*exp(-i_*kz*(z-zs));
  }
  else if (layi == 0) {
    kz = calcKz(0,true);
    aSi.t[1] = S.t[0][1]*(aS.t[1]+a.t[1])*exp(-i_*kz*(z-Zi.t[0]));
    aSi.t[0] = 0.;
  }
  else if (layi == nbrMediums-1) {
    S = calcSMatrixU(pol,zs);
    kz = calcKz(nbrMediums-1,true);
    aSi.t[0] = S.t[1][0]*(aS.t[0]+a.t[0])*exp(i_*kz*(z-Zi.t[nbrMediums-2]));
    aSi.t[1] = 0.;
  }
  else if (z<zs) {
    aSi.t[0] = (aS.t[0]-S.t[1][1]*(aS.t[1]+a.t[1]))/S.t[1][0];
    aSi.t[1] = (S.t[0][0]*aS.t[0]+(S.t[0][1]*S.t[1][0]-S.t[0][0]*S.t[1][1])*(aS.t[1]+a.t[1]))/S.t[1][0];
  }
  else if (z>zs) {
    aSi.t[0] = ((S.t[0][1]*S.t[1][0]-S.t[0][0]*S.t[1][1])*(aS.t[0]+a.t[0])+S.t[1][1]*aS.t[1])/S.t[0][1];
    aSi.t[1] = (aS.t[1]-S.t[0][0]*(aS.t[0]+a.t[0]))/S.t[0][1];
  }
  aSi = exp(i_*(kx*pos.t[0]+ky*pos.t[1]))*aSi;

  kz = calcKz(layi,true);
  Complex bP = acos(kz/ki.t[layi]), bM = acos(-kz/ki.t[layi]);
  Vector3d EP(0.),EM(0.);

  EP.t[2-pol] = aSi.t[0];
  EM.t[2-pol] = aSi.t[1];
  EP = toCartesianCoordinates(EP,bP,alpha0);
  EM = toCartesianCoordinates(EM,bM,alpha0);

  return EM+EP;
}

Vector Multilayer::calcSourceAmplIn(const double &z, const int &pol) {
  Vector aPM(2);
  Matrix Sl,Su;
  Sl = calcSMatrixL(pol,z);
  Su = calcSMatrixU(pol,z);
  Complex cst = 1.-Su.t[0][0]*Sl.t[1][1];
  aPM.t[0] = Sl.t[1][1]*(Su.t[0][0]*a.t[0]+a.t[1])/cst;
  aPM.t[1] = Su.t[0][0]*(Sl.t[1][1]*a.t[1]+a.t[0])/cst;
  return aPM;
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

void Multilayer::setLambda(const double &lambd) {
  lambda = lambd;
  int i;
  for (i=0 ; i<nbrMediums ; i++) {
    ni.t[i] = ri[i].calcIndex(lambda);
    ki.t[i] = 2.*piDbl*ni.t[i]/lambda;
  }
  if (beta0<piDbl/2.) {
    kx = ki.t[0]*cos(alpha0)*sin(beta0);
    ky = ki.t[0]*sin(alpha0)*sin(beta0);
  }
  else {
    kx = ki.t[nbrMediums-1]*cos(alpha0)*sin(beta0);
    ky = ki.t[nbrMediums-1]*sin(alpha0)*sin(beta0);
  }
  return;
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
  if (lay == 0) return calcSMatrix(pol);
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
  if (lay == nbrMediums-1) return calcSMatrix(pol);
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
  if (lay1 == 0) return calcSMatrixL(pol,z2);
  if (lay2 == nbrMediums-1) return calcSMatrixU(pol,z1);
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




