#ifndef MULTIPLESCATTERING_H_INCLUDED
#define MULTIPLESCATTERING_H_INCLUDED

#include "AxisSymParticle.h"
#include "Translation.h"
#include "IncidentRadiation.h"
#include "Reflection.h"
#include "Multilayer.h"

class MultipleScattering
{
public:
  MultipleScattering(const string &fname);
  void calcAmni();
  void calcTAmni();
  void calcAmniElectron();
  void calcTAmniElectron();
  void calcPmn();
  Vector calcPmni0();
  void solve();
  void gmres(Vector &X, const double &tol);
  Matrix calcDirect();
  Vector calcDirectGMRES(const Vector &B);
  Vector calcFFT1D(const Vector &B);
  void constructMatrix1D();
  Vector calcFFT2D(const Vector &B);
  void constructMatrix2D();

  void setN(const int &n);
  void setLambda(const double &l);
  void setElectronEnergy(const double &ee);

  double calcCsca(const bool &normalize, const double &nor);
  RVector calcCabs(const bool &normalize, const double &nor);
  RVector calcCext(const bool &normalize, const double &nor);
  RVector calcCextOT(const bool &normalize, const double &nor);
  RVector calcScaML(const bool &normalize, const double &nor);
  double calcCscaN(const bool &normalize, const double &nor, const int &nn, const int &pol);
  RVector calcCextN(const bool &normalize, const double &nor, const int &nn, const int &pol);
  RVector calcCabsN(const bool &normalize, const double &nor, const int &nn, const int &pol);
  Vector calcComplexCext(const bool &normalize, const double &nor);
  Vector calcComplexCextN(const bool &normalize, const double &nor, const int &nn, const int &pol);
  Vector calcStokesParametersInc();
  Vector calcStokesParametersExt();
  Matrix calcScatteringMatrix(const double &ctp, const double &ctm, const double &alphap, const double &alpham);
  Matrix calcExtinctionMatrix();

  Vector3d calcEsca(const double &xx, const double &yy, const double &zz);
  Vector3d calcHsca(const double &xx, const double &yy, const double &zz);
  Vector3d calcEscaInf(const double &theta, const double &phi);
  Vector3d calcHscaInf(const double &theta, const double &phi);
  Vector3d calcEint(const int &part, const double &xx, const double &yy, const double &zz);
  Vector3d calcHint(const int &part, const double &xx, const double &yy, const double &zz);
  Vector3d calcEinc(const double &xx, const double &yy, const double &zz);
  Vector3d calcHinc(const double &xx, const double &yy, const double &zz);

  double calcEEL();
  double calcCL();
  double calcEEL(const int &n, const int &m, const int &mode);

  AxisSymParticle *particles;
  IncidentRadiation incidence;
  Translation trans;
  Reflection refl;
  Reflection1 refl1;
  Multilayer1 ml;
  int Np,Npx,Npy,Npz,N,NM,solver,Ng,methodTMatrix,periodic,nbrDim,nMatFFT,nMatFFTx,nMatFFTy,lay,nMatrices;
  double *x,*y,*z,xc,yc,zc,tol,lambda,dx,dy,dz,z0;
  Vector Amni,AmnRi,AmnTi,Pmni,Cmni,Pmn;
  Matrix *matricesFFT1D;
  RMatrix positions;
  Complex nm1,nm2;
};

#endif // MULTIPLESCATTERING_H_INCLUDED
