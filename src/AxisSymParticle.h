#ifndef AXISSYMPARTICLE_H_INCLUDED
#define AXISSYMPARTICLE_H_INCLUDED

#include "RefractiveIndex.h"
#include "RealContainers.h"
#include "Matrix.h"
#include "BesselFunctions.h"
#include "LegendreFunctions.h"
#include "Integration.h"
#include "VSWF.h"

class AxisSymParticle
{
public:
    AxisSymParticle();
    void setParticleProperties(const int &typ,const int &mat,const int &matExt,const double &lambd,const RVector &rr);
    void setLambda(const double &lambd);
    double calcR(const double &cost);
    double calcDR(const double &cost);
    double lambda;
    RVector geom;
    int type,material,materialExt;
    Complex n1,k1,n2,k2;

    void setTMatrixParameters(const int &n,const int &ng,const int &meth);
    void calcTSphere();
    void calcTEBCM();
    void calcTLS();
    void calcTBC();
    void calcTLScylinder();
    void calcSources();
    void calcTDS();
    Matrix getT();
    Matrix getR();
    Vector multT(const Vector &v);
    int N,Ng,method,nbr;
    Matrix *T,*R,*Refl;
    Vector sources;
};

#endif // AXISSYMPARTICLE_H_INCLUDED
