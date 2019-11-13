#ifndef REFLECTION_H_INCLUDED
#define REFLECTION_H_INCLUDED

#include "Matrix.h"
#include "BesselFunctions.h"
#include "LegendreFunctions.h"
#include "Integration.h"
#include "Multilayer.h"
#include "Multilayer1.h"
#include "ClebshGordan.h"
#include "WignerFunctions.h"

class Reflection
{
public:
  Reflection();
  void setParameters(const double &xx, const double &yy, const double &zzi, const double &zzj, const bool &self);
  //void setVSWFParameters(const double &xx, const double &yy, const double &zzi, const double zzj);
  void setConstants(const double &lambd, const int &n, const int &ng);
  //void calcVSWFsIntegral();
  Matrix calcR();
  Matrix calcRSelf();
  void calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj);

  int N,NM,Ng,lay,layI,layJ,nbrLay;
  double dx,dy,zi,zj,xMin,xMax,yMax,rauij,lambda;
  Complex eiphi;
  Matrix R,Ri,legPi,legTau,gaussLeg,bessj,Mmn,Nmn;
  Vector kappa,dkappa,kz;
  Multilayer1 ml;
  bool Self;
  LegendreFunctionsC leg;
  Bessel bes;
};

class Reflection1
{
public:
  Reflection1();
  void setParameters(const double &xx, const double &yy, const double &zzi, const double zzj, const bool &self);
  void setPositions(const RMatrix &pos);
  //void setVSWFParameters(const double &xx, const double &yy, const double &zzi, const double zzj);
  void setConstants(const double &lambd, const int &n, const int &ng);
  //void calcVSWFsIntegral();
  void calcR();
  void calcRSelf();
  int nbrQelmts();
  Complex getQelmt(const int &u, const int &v, const int &w, const int &pt, const int &pm);
  int getElmt(const int &u, const int &v, const int &w);
  int getRefl(const int &i, const int &j);
  int getNbrRefl();
  void calcQintegrals();
  void calcQintegrals(const int &i, const int &j);
  void calcQintegralsSelf();
  void calcQintegralsSelf(const int &i, const int &j);
  void calcRMatrices(const int &i, const int &j);
  void calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj);

  int N,NM,Ng,lay,layI,layJ,nbrLay,nbrQ,nbrReflections;
  double dx,dy,zi,zj,xMin,xMax,yMax,rauij,lambda;
  Complex eiphi,eiphii;
  Matrix R,Ri,legPi,legTau,gaussLeg,bessj,duvw,Qintegrals,Mmn,Nmn;
  RMatrix positions,gaussLegAlpha,gaussLegBeta1,gaussLegBeta2;
  Vector kappa,dkappa,kz;
  Multilayer1 ml;
  bool Self;
  LegendreFunctionsC leg;
  Bessel bes;
  Vector Q1,Q2,Q3,Q4;
  WignerFunctionsC wign;
  ClebshGordan cg;
};

#endif // REFLECTION_H_INCLUDED
