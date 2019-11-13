#ifndef INCIDENTRADIATION_H_INCLUDED
#define INCIDENTRADIATION_H_INCLUDED

#include "Vector.h"
#include "LegendreFunctions.h"
#include "BesselFunctions.h"
#include "SpecialFunctions.h"
#include "Multilayer1.h"
#include <fstream>

class IncidentRadiation
{
public:
  IncidentRadiation();
  IncidentRadiation(const int &inc, const double &lambd);
  void setParameters(const int &inc, const double &lambd);
  void setPlaneWaveCoeffs(const int &n, const Vector3d &pos);
  void setElectron(const int &n, const double &E0, const double &E, const double &x, const double &y);

  int N,type,incidence;
  Multilayer1 ml;
  // plane wave
  double beta0,alpha0,sinb0;
  Complex Eb, Ea, cosb0;
  // electron
  double E0,E,v0,omega,gamma;
  double x0,y0,z0,b,phi0;
  double lambda;
  Vector aincs,aincp,aouts,aoutp,ap,as;

  Vector Amn,AmnR,AmnT;
  string fileName;
};

#endif // INCIDENTRADIATION_H_INCLUDED
