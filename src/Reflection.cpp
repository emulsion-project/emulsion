#include "Reflection.h"
#include <iostream>
using namespace std;

Reflection::Reflection() {
  ml.setParameters("Multilayer/multilayerConfig.txt",488e-9);
  N = 0;
  Ng = 0;
  lambda = 0.;
}

void Reflection::setConstants(const double &lambd, const int &n, const int &ng) {
  if ((N != n) || (Ng != ng) || (lambd !=lambda)) {
    lambda = lambd;
    Ng = ng;
    N = n;
    NM = (N+1)*(N+1)-1;
    ml.setLambda(lambd);
    nbrLay = ml.nbrMediums;
    xMin = 0.;
    for (int i=0 ; i<nbrLay ; i++) if (xMin < real(ml.ki.t[i])) xMin = real(ml.ki.t[i]);
    xMin += 2.*piDbl/lambd;
    xMax = xMin*200.;
    yMax = xMin/10.;

    RMatrix gaussLeg1 = gaussLegQuad(Ng,0.,piDbl);
    RMatrix gaussLeg2 = gaussLegQuad(Ng,xMin,xMax);
    gaussLeg = Matrix(2*nbrLay*Ng,2);

    legPi = Matrix(2*nbrLay*Ng,NM);
    legTau = Matrix(2*nbrLay*Ng,NM);
    kappa = Vector(2*Ng);
    dkappa = Vector(2*Ng);
    kz = Vector(2*nbrLay*Ng);
    double cstn;

    for (int l=0 ; l<nbrLay ; l++) {

      for (int i=0 ; i<Ng ; i++) {

        gaussLeg.t[i][0] = gaussLeg1.t[i][0];
        gaussLeg.t[i][1] = gaussLeg1.t[i][1];
        gaussLeg.t[i+Ng][0] = gaussLeg2.t[i][0];//
        gaussLeg.t[i+Ng][1] = gaussLeg2.t[i][1];//

        kappa.t[i] = xMin*0.5*(1.-cos(gaussLeg1.t[i][0])) - i_*yMax*sin(gaussLeg1.t[i][0]);
        dkappa.t[i] = xMin*0.5*sin(gaussLeg1.t[i][0]) - i_*yMax*cos(gaussLeg1.t[i][0]);
        kappa.t[i+Ng] = gaussLeg2.t[i][0];//
        dkappa.t[i+Ng] = 1.;

        ml.setKxy(kappa.t[i]);
        kz.t[i+2*l*Ng] = ml.calcKz(l,true);
        leg.assign(N,kz.t[i+2*l*Ng]/ml.ki.t[l]);
        for (int nn=1 ; nn<=N ; nn++) {
          cstn = 1./sqrt(double(2*nn*(nn+1)));
          for (int m=-nn ; m<=nn ; m++) {
            legPi.t[i+2*l*Ng][ii(nn,m)] = cstn*double(m)*leg.Pimn(nn,m);
            legTau.t[i+2*l*Ng][ii(nn,m)] = cstn*leg.Taumn(nn,m);
          }
        }

        ml.setKxy(kappa.t[i+Ng]);
        kz.t[i+Ng+2*l*Ng] = ml.calcKz(l,true);
        leg.assign(N,kz.t[i+Ng+2*l*Ng]/ml.ki.t[l]);
        for (int nn=1 ; nn<=N ; nn++) {
          cstn = 1./sqrt(double(2*nn*(nn+1)));
          for (int m=-nn ; m<=nn ; m++) {
            legPi.t[i+Ng+2*l*Ng][ii(nn,m)] = cstn*double(m)*leg.Pimn(nn,m);
            legTau.t[i+Ng+2*l*Ng][ii(nn,m)] = cstn*leg.Taumn(nn,m);
          }
        }
      }
    }

  }
  return;
}

void Reflection::setParameters(const double &xx, const double &yy, const double &zzi, const double &zzj, const bool &self) {
  lay = ml.getMedium(zzi);
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  dx = xx;
  dy = yy;
  zi = zzi;
  zj = zzj;

  if (!self) {
    bessj = Matrix(2*N+1,2*Ng);
    rauij = sqrt(dx*dx+dy*dy);
    eiphi = (dx+i_*dy)/rauij;
    if (rauij == 0.) eiphi = 1.;
    for (int kk=0 ; kk<Ng ; kk++) {
      bes.assign(2*N,rauij*kappa.t[kk]);
      for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
    }
    for (int kk=Ng ; kk<2*Ng ; kk++) {
      bes.assign(2*N,rauij*kappa.t[kk]);
      for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
    }
  }
  if (self) R = calcRSelf();
  else R = calcR();
  return;
}

Matrix Reflection::calcRSelf() {
  Matrix  RR(2*NM,2*NM,0.);
  int line=0,col=0;
  Complex cstMult;
  Complex mult,PP,TT,TP,PT;
  Complex rppp,rppm,rptp,rptm,rspp,rspm,rstp,rstm;
  Vector app,asp,apt,ast,a(2,1.);
  for (int k=0 ; k<2*Ng ; k++) {

    ml.setKxy(kappa.t[k]);
    cstMult = gaussLeg.t[k][1]*4.*kappa.t[k]*dkappa.t[k]/ml.ki.t[lay]/kz.t[k+2*lay*Ng];
    line = 0;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        a.t[1] = powm1(n-m);
        app = ml.sourceInt(1,a,zj,zj);
        asp = ml.sourceInt(0,a,zj,zj);
        a.t[1] = -powm1(n-m);
        apt = ml.sourceInt(1,a,zj,zj);
        ast = ml.sourceInt(0,a,zj,zj);



        col=0;
        for (int l=1 ; l<=N ; l++) {
          for (int kk=-l ; kk<=l ; kk++) {
            if (m == kk) {
              mult = powi(l-n)*cstMult;
              rppp = app.t[0] + powm1(l-kk)*app.t[1];
              rppm = app.t[0] - powm1(l-kk)*app.t[1];
              rptp = apt.t[0] + powm1(l-kk)*apt.t[1];
              rptm = apt.t[0] - powm1(l-kk)*apt.t[1];
              rspp = asp.t[0] + powm1(l-kk)*asp.t[1];
              rspm = asp.t[0] - powm1(l-kk)*asp.t[1];
              rstp = ast.t[0] + powm1(l-kk)*ast.t[1];
              rstm = ast.t[0] - powm1(l-kk)*ast.t[1];

              PP = legPi.t[k+2*lay*Ng][ii(n,m)]*legPi.t[k+2*lay*Ng][ii(l,m)];
              TT = legTau.t[k+2*lay*Ng][ii(n,m)]*legTau.t[k+2*lay*Ng][ii(l,m)];
              PT = legPi.t[k+2*lay*Ng][ii(n,m)]*legTau.t[k+2*lay*Ng][ii(l,m)];
              TP = legTau.t[k+2*lay*Ng][ii(n,m)]*legPi.t[k+2*lay*Ng][ii(l,m)];

              RR.t[line][col] += mult*(PP*rppp+TT*rstm);
              RR.t[line][col+NM] += mult*(PT*rppm+TP*rstp);
              RR.t[line+NM][col] += mult*(TP*rptp+PT*rspm);
              RR.t[line+NM][col+NM] += mult*(TT*rptm+PP*rspp);
            }
            col++;
          }
        }
        line++;
      }
    }
  }
  return RR;
}

Matrix Reflection::calcR() {
  Matrix  RR(2*NM,2*NM,0.);
  int line=0,col=0;
  Vector app,asp,apt,ast,a(2,1.);
  Vector vv(3*Ng);
  Complex mult,rpp,rpm,rsp,rsm,PP,TT,TP,PT,cstMult;
  Complex rppp,rppm,rptp,rptm,rspp,rspm,rstp,rstm;
  for (int k=0 ; k<2*Ng ; k++) {
    ml.setKxy(kappa.t[k]);
    cstMult = gaussLeg.t[k][1]*4.*kappa.t[k]*dkappa.t[k]/ml.ki.t[layI]/kz.t[k+2*layI*Ng];
    line = 0;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        a.t[1] = powm1(n-m);
        app = ml.sourceInt(1,a,zj,zi);
        asp = ml.sourceInt(0,a,zj,zi);
        a.t[1] = -powm1(n-m);
        apt = ml.sourceInt(1,a,zj,zi);
        ast = ml.sourceInt(0,a,zj,zi);

        //app.t[1] = asp.t[1] = powm1(n-m)*exp(i_*kz.t[k+2*lay*Ng]*(zi-zj));
        //apt.t[1] = ast.t[1] = -powm1(n-m)*exp(i_*kz.t[k+2*lay*Ng]*(zi-zj));
        //app.t[0] = asp.t[0] = apt.t[0] = ast.t[0] = 0.;
        col = 0;
        for (int l=1 ; l<=N ; l++) {
          for (int kk=-l ; kk<=l ; kk++) {
              mult = powi(l-n+abs(m-kk))*pow(eiphi,double(m-kk))*bessj.t[int(abs(m-kk))][k]*cstMult;
              rppp = app.t[0] + powm1(l-kk)*app.t[1];
              rppm = app.t[0] - powm1(l-kk)*app.t[1];
              rptp = apt.t[0] + powm1(l-kk)*apt.t[1];
              rptm = apt.t[0] - powm1(l-kk)*apt.t[1];
              rspp = asp.t[0] + powm1(l-kk)*asp.t[1];
              rspm = asp.t[0] - powm1(l-kk)*asp.t[1];
              rstp = ast.t[0] + powm1(l-kk)*ast.t[1];
              rstm = ast.t[0] - powm1(l-kk)*ast.t[1];

              PP = legPi.t[k+2*layJ*Ng][ii(n,m)]*legPi.t[k+2*layI*Ng][ii(l,kk)];
              TT = legTau.t[k+2*layJ*Ng][ii(n,m)]*legTau.t[k+2*layI*Ng][ii(l,kk)];
              PT = legPi.t[k+2*layJ*Ng][ii(n,m)]*legTau.t[k+2*layI*Ng][ii(l,kk)];
              TP = legTau.t[k+2*layJ*Ng][ii(n,m)]*legPi.t[k+2*layI*Ng][ii(l,kk)];

              RR.t[line][col] += mult*(PP*rppp+TT*rstm);
              RR.t[line][col+NM] += mult*(PT*rppm+TP*rstp);
              RR.t[line+NM][col] += mult*(TP*rptp+PT*rspm);
              RR.t[line+NM][col+NM] += mult*(TT*rptm+PP*rspp);
            col++;
          }
        }
        line++;
      }
    }
  }
  return RR;
}

/*void Reflection::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector app,asp,apt,ast,a(2,1.);
  Complex rpp,rpt,rsp,rst;
  Complex cost,mult,multExp,multExp1,multExp2;
  double alpha = atan2(dy,dx);
  int elmtJ,elmtI;
  Complex cosb;
  rauij = sqrt(dx*dx+dy*dy);
  eiphi = (dx+i_*dy)/rauij;
  if (rauij==0.) eiphi=1.;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  //ml.setLambda(lambda);
  RMatrix legAlpha = gaussLegQuad(Ng/2,0.,2.*piDbl);
  for (int k=5 ; k<6 ; k++) {
    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    bes.assign(N+1,rauij*kappa.t[k]);
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[elmtI]/ml.ki.t[layI];

    for (int l=5 ; l<6 ; l++) {
      int i=0;
      alpha = legAlpha.t[l][0];
      for (int n=1 ; n<=N ; n++) {
        for (int m=-n ; m<=n ; m++) {

          a.t[0] = 1.;
          a.t[1] = powm1(n-m);
          a.t[1] = 0.;
          app = ml.sourceInt(1,a,zzj,zzi);
          asp = ml.sourceInt(0,a,zzj,zzi);
          if (layI==layJ) {
            app = app + exp(i_*ml.calcKz(layI,true)*fabs(zzj-zzi))*a;
            asp = asp + exp(i_*ml.calcKz(layI,true)*fabs(zzj-zzi))*a;
          }
          a.t[0] = 1.;
          a.t[1] = -powm1(n-m);
          a.t[1] = 0.;

          apt = ml.sourceInt(1,a,zzj,zzi);
          ast = ml.sourceInt(0,a,zzj,zzi);
          if (layI==layJ) {
            apt = apt + exp(i_*ml.calcKz(layI,true)*fabs(zzj-zzi))*a;
            ast = ast + exp(i_*ml.calcKz(layI,true)*fabs(zzj-zzi))*a;
          }

          Mmn.t[i][0] += mult*legAlpha.t[l][1]*((app.t[0]-app.t[1])*cosb*cos(alpha) - (apt.t[0]-apt.t[1])*sin(alpha));
          Mmn.t[i][1] += mult*legAlpha.t[l][1]*((app.t[0]-app.t[1])*cosb*sin(alpha) + (apt.t[0]+apt.t[1])*cos(alpha));
          Mmn.t[i][2] += mult*legAlpha.t[l][1]*((app.t[0]-app.t[1])*sqrt(1.-cosb*cosb));

          Nmn.t[i][0] += mult*legAlpha.t[l][1]*((apt.t[0]+apt.t[1])*cosb*cos(alpha) - (app.t[0]-app.t[1])*sin(alpha));
          Nmn.t[i][1] += mult*legAlpha.t[l][1]*((apt.t[0]+apt.t[1])*cosb*sin(alpha) + (app.t[0]-app.t[1])*cos(alpha));
          Nmn.t[i][2] += mult*legAlpha.t[l][1]*((apt.t[0]-apt.t[1])*sqrt(1.-cosb*cosb));
          i++;
        }
      }
    }
  }
  return;
}*/

void Reflection::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector app,asp,apt,ast,a(2,1.);
  Complex rpp,rpt,rsp,rst;
  Complex cost,mult,multExp,multExp1,multExp2;
  double alpha = atan2(dy,dx);
  int elmtJ,elmtI;
  Complex cosb;
  rauij = sqrt(dx*dx+dy*dy);
  eiphi = (dx+i_*dy)/rauij;
  if (rauij==0.) eiphi=1.;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  //ml.setLambda(lambda);
  for (int k=0 ; k<2*Ng ; k++) {

    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    bes.assign(N+1,rauij*kappa.t[k]);
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[elmtI]/ml.ki.t[layI];
    int i=0;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        a.t[0] = 1.;
        a.t[1] = powm1(n-m);
        app = ml.sourceInt(1,a,zzj,zzi);
        asp = ml.sourceInt(0,a,zzj,zzi);
        a.t[0] = 1.;
        a.t[1] = -powm1(n-m);
        apt = ml.sourceInt(1,a,zzj,zzi);
        ast = ml.sourceInt(0,a,zzj,zzi);

        multExp = mult*bes.J.t[int(abs(m))]*powi(-n+abs(m)-1)*pow(eiphi,double(m));
        multExp1 = mult*bes.J.t[int(abs(m+1))]*powi(-n+abs(m+1)-1)*pow(eiphi,double(m+1));
        multExp2 = mult*bes.J.t[int(abs(m-1))]*powi(-n+abs(m-1)-1)*pow(eiphi,double(m-1));

        Mmn.t[i][0] += 0.5*((multExp1+multExp2)*legPi.t[elmtJ][i]*(app.t[0]-app.t[1])*cosb - (multExp1-multExp2)*legTau.t[elmtJ][i]*(ast.t[0]+ast.t[1]));
        Mmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*legPi.t[elmtJ][i]*(app.t[0]-app.t[1])*cosb + (multExp1+multExp2)*i_*legTau.t[elmtJ][i]*(ast.t[0]+ast.t[1]));
        Mmn.t[i][2] += multExp*((app.t[0]+app.t[1])*legPi.t[elmtJ][i]*sqrt(1.-cosb*cosb));

        Nmn.t[i][0] += 0.5*((multExp1+multExp2)*legTau.t[elmtJ][i]*(apt.t[0]-apt.t[1])*cosb - (multExp1-multExp2)*legPi.t[elmtJ][i]*(asp.t[0]+asp.t[1]));
        Nmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*legTau.t[elmtJ][i]*(apt.t[0]-apt.t[1])*cosb + (multExp1+multExp2)*i_*legPi.t[elmtJ][i]*(asp.t[0]+asp.t[1]));
        Nmn.t[i][2] += multExp*(legTau.t[elmtJ][i]*(apt.t[0]+apt.t[1])*sqrt(1.-cosb*cosb));
        i++;
      }
    }
  }
  return;
}

/*void Reflection::calcVSWFsIntegral() {
  for (int i=0 ; i<3 ; i++) Mmn[i].t[i] = Nmn[i].t[i] = 0.;
  Vector as(2),ap(2);
  Vector3d pos;
  pos.t[0] = 0.;
  pos.t[1] = 0.;
  pos.t[2] = zi;
  Complex fact,cost,mult,tp,ts,q,multExp;
  double alpha;

  for (int ka=0 ; ka<Ng ; ka++) {
    alpha = gaussLegP.t[ka][0];

    for (int kb=0 ; kb<Ng ; kb++) {
      Complex jfactor = ml.ni.t[layI]*ml.ni.t[layI]*gaussLegT.t[kb][0]/(ml.ni.t[layJ]*ml.ni.t[layJ]*yInt[kb]);
      ap = ml.setSourceIntAmpl(1,1.,1.,acos(yInt[kb]),0.,zj,pos);
      as = ml.setSourceIntAmpl(0,1.,1.,acos(yInt[kb]),0.,zj,pos);
      //ap.t[0] = as.t[0] = 0.;
      //as.t[1] = ap.t[1] = 1.;

      cost = gaussLegT.t[kb][0];
      mult = -exp(i_*ml.ki.t[0]*((dx*cos(alpha)+dy*sin(alpha))*sqrt(1.-cost*cost)+zi*cost));
      mult *= gaussLegP.t[ka][1]*gaussLegT.t[kb][1]*jfactor;
      //cout << jfactor << endl;
      int i=0;
      for (int n=1 ; n<=N ; n++) {
        for (int m=-n ; m<=n ; m++) {
          multExp = exp(i_*double(m)*alpha)/(2.*piDbl*powi(n+1)*sqrt(double(2*n*(n+1))));
          Mmn[i].t[0] += multExp*mult*(double(m)*legS[kb].Pimn(n,m)*ap.t[1]*(-cost)*cos(alpha)-i_*legS[kb].Taumn(n,m)*as.t[1]*sin(alpha));
          Mmn[i].t[1] += multExp*mult*(double(m)*legS[kb].Pimn(n,m)*ap.t[1]*(-cost)*sin(alpha)+i_*legS[kb].Taumn(n,m)*as.t[1]*cos(alpha));
          Mmn[i].t[2] -= multExp*mult*double(m)*legS[kb].Pimn(n,m)*ap.t[1]*sqrt(1.-cost*cost);

          Mmn[i].t[0] += multExp*mult*powm1(m-n+1)*(double(m)*legS[kb].Pimn(n,m)*ap.t[0]*(-cost)*cos(alpha)-i_*legS[kb].Taumn(n,m)*as.t[0]*sin(alpha));
          Mmn[i].t[1] += multExp*mult*powm1(m-n+1)*(double(m)*legS[kb].Pimn(n,m)*ap.t[0]*(-cost)*sin(alpha)+i_*legS[kb].Taumn(n,m)*as.t[0]*cos(alpha));
          Mmn[i].t[2] -= multExp*mult*powm1(m-n+1)*double(m)*legS[kb].Pimn(n,m)*ap.t[0]*sqrt(1.-cost*cost);

          Nmn[i].t[0] += multExp*mult*(legS[kb].Taumn(n,m)*ap.t[1]*(-cost)*cos(alpha)-i_*double(m)*legS[kb].Pimn(n,m)*as.t[1]*sin(alpha));
          Nmn[i].t[1] += multExp*mult*(legS[kb].Taumn(n,m)*ap.t[1]*(-cost)*sin(alpha)+i_*double(m)*legS[kb].Pimn(n,m)*as.t[1]*cos(alpha));
          Nmn[i].t[2] -= multExp*mult*legS[kb].Taumn(n,m)*ap.t[1]*sqrt(1.-cost*cost);

          Nmn[i].t[0] += multExp*mult*powm1(m-n+1)*(legS[kb].Taumn(n,m)*ap.t[0]*(-cost)*cos(alpha)-i_*double(m)*legS[kb].Pimn(n,m)*as.t[0]*sin(alpha));
          Nmn[i].t[1] += multExp*mult*powm1(m-n+1)*(legS[kb].Taumn(n,m)*ap.t[0]*(-cost)*sin(alpha)+i_*double(m)*legS[kb].Pimn(n,m)*as.t[0]*cos(alpha));
          Nmn[i].t[2] -= multExp*mult*powm1(m-n+1)*legS[kb].Taumn(n,m)*ap.t[0]*sqrt(1.-cost*cost);
          i++;
        }
      }
    }
  }
  return;
}*/

Reflection1::Reflection1() {
  ml.setParameters("Multilayer/multilayerConfig.txt",488e-9);
  N = 0;
  Ng = 0;
  lambda = 0.;
}

int Reflection1::nbrQelmts() {
  int i=0;
  for (int u=-2 ; u<=2 ; u+=2) {
    for (int w=abs(u) ; w<=2*N ; w++) {
      for (int v=0 ; v<=w ; v++) {
        i++;
      }
    }
  }
  return i;
}



int Reflection1::getElmt(const int &u, const int &v, const int &w) {
  int i=0,elmt;
  for (int U=-2 ; U<=2 ; U+=2) {
    for (int W=abs(U) ; W<=2*N ; W++) {
      for (int V=0 ; V<=W ; V++) {
        if(U==u && V==v && W==w) return i;
        i++;
      }
    }
  }
  return i-1;
}

void Reflection1::setConstants(const double &lambd, const int &n, const int &ng) {
  if (N != n) cg.assign(n);
  if ((N != n) || (Ng != ng) || (lambd !=lambda)) {
    Ng = ng;
    N = n;
    NM = (N+1)*(N+1)-1;
    ml.setLambda(lambd);
    lambda = lambd;
    nbrLay = ml.nbrMediums;
    nbrQ = nbrQelmts();
    xMin = 0.;
    for (int i=0 ; i<nbrLay ; i++) if (xMin < real(ml.ki.t[i])) xMin = real(ml.ki.t[i]);
    xMin += 2.*piDbl/lambd;
    xMax = xMin*100.;
    yMax = xMin/10.;

    RMatrix gaussLeg1 = gaussLegQuad(Ng,0.,piDbl);
    gaussLegAlpha = gaussLegQuad(Ng,0.,2.*piDbl);
    RMatrix gaussLeg2 = gaussLegQuad(Ng,xMin,xMax);
    gaussLeg = Matrix(2*nbrLay*Ng,2);
    duvw = Matrix(2*nbrLay*Ng,nbrQ);
    legPi = Matrix(2*nbrLay*Ng,NM);
    legTau = Matrix(2*nbrLay*Ng,NM);
    kappa = Vector(2*Ng);
    dkappa = Vector(2*Ng);
    kz = Vector(2*nbrLay*Ng);
    int elmt;
    double cst,cstn;
    for (int l=0 ; l<nbrLay ; l++) {
      for (int i=0 ; i<Ng ; i++) {
        gaussLeg.t[i][0] = gaussLeg1.t[i][0];
        gaussLeg.t[i][1] = gaussLeg1.t[i][1];
        gaussLeg.t[i+Ng][0] = gaussLeg2.t[i][0];//
        gaussLeg.t[i+Ng][1] = gaussLeg2.t[i][1];//

        kappa.t[i] = xMin*0.5*(1.-cos(gaussLeg1.t[i][0])) - i_*yMax*sin(gaussLeg1.t[i][0]);
        dkappa.t[i] = xMin*0.5*sin(gaussLeg1.t[i][0]) - i_*yMax*cos(gaussLeg1.t[i][0]);
        kappa.t[i+Ng] = gaussLeg2.t[i][0];//
        dkappa.t[i+Ng] = 1.;

        ml.setKxy(kappa.t[i]);
        kz.t[i+2*l*Ng] = ml.calcKz(l,true);
        wign.assign(2*N,kz.t[i+2*l*Ng]/ml.ki.t[l]);
        elmt=0;
        for (int u=-2 ; u<=2 ; u+=2) {
          for (int w=abs(u) ; w<=2*N ; w++) {
            for (int v=0 ; v<=w ; v++) {
              duvw.t[i+2*l*Ng][elmt] = wign.dlmn(u,v,w);
              elmt++;
            }
          }
        }

        leg.assign(N,kz.t[i+2*l*Ng]/ml.ki.t[l]);
        for (int nn=1 ; nn<=N ; nn++) {
          cstn = 1./sqrt(double(2*nn*(nn+1)));
          for (int m=-nn ; m<=nn ; m++) {
            legPi.t[i+2*l*Ng][ii(nn,m)] = cstn*double(m)*leg.Pimn(nn,m);
            legTau.t[i+2*l*Ng][ii(nn,m)] = cstn*leg.Taumn(nn,m);
          }
        }

        ml.setKxy(kappa.t[i+Ng]);
        kz.t[i+Ng+2*l*Ng] = ml.calcKz(l,true);
        wign.assign(2*N,kz.t[i+Ng+2*l*Ng]/ml.ki.t[l]);
        elmt=0;
        for (int u=-2 ; u<=2 ; u+=2) {
          for (int w=abs(u) ; w<=2*N ; w++) {
            for (int v=0 ; v<=w ; v++) {
              duvw.t[i+Ng+2*l*Ng][elmt] = wign.dlmn(u,v,w);
              elmt++;
            }
          }
        }
        leg.assign(N,kz.t[i+Ng+2*l*Ng]/ml.ki.t[l]);
        for (int nn=1 ; nn<=N ; nn++) {
          cstn = 1./sqrt(double(2*nn*(nn+1)));
          for (int m=-nn ; m<=nn ; m++) {
            legPi.t[i+Ng+2*l*Ng][ii(nn,m)] = cstn*double(m)*leg.Pimn(nn,m);
            legTau.t[i+Ng+2*l*Ng][ii(nn,m)] = cstn*leg.Taumn(nn,m);
          }
        }
      }
    }

  }
  return;
}

void Reflection1::setParameters(const double &xx, const double &yy, const double &zzi, const double zzj, const bool &self) {
  lay = ml.getMedium(zi);
  layI = ml.getMedium(zi);
  layJ = ml.getMedium(zj);
  dx = -xx;
  dy = -yy;
  zi = zzi;
  zj = zzj;

  if (!self) {
    bessj = Matrix(2*N+1,2*Ng);
    rauij = sqrt(dx*dx+dy*dy);
    eiphi = (dx+i_*dy)/rauij;
    if (rauij == 0.) eiphi = 1.;
    for (int kk=0 ; kk<Ng ; kk++) {
      bes.assign(2*N,rauij*kappa.t[kk]);
      for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
    }
    for (int kk=Ng ; kk<2*Ng ; kk++) {
      bes.assign(2*N,rauij*kappa.t[kk]);
      for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
    }
  }

  if (self) {
    calcQintegralsSelf();
    calcRSelf();
  }
  else {
    calcQintegrals();
    calcR();
  }
  return;
}

int Reflection1::getNbrRefl() {
  int nbr=positions.n*4;
  for (int i=0 ; i<positions.n ; i++) {
    for (int j=i+1 ; j<positions.n ; j++) {
      nbr+=4;
    }
  }
  return nbr;
}

int Reflection1::getRefl(const int &i, const int &j) {
  if (i==j) return i*4;
  int ii=i,jj=j,nbr=positions.n*4,refl;
  if (ii>jj) swap(ii,jj);
  for (int a=0 ; a<positions.n ; a++) {
    for (int b=a+1 ; b<positions.n ; b++) {
      if (a==ii && b==jj) refl=nbr;
      nbr+=4;
    }
  }
  return refl;
}

void Reflection1::setPositions(const RMatrix &pos) {
  positions = pos;
  nbrReflections = getNbrRefl();
  Qintegrals = Matrix(nbrQ,nbrReflections);
  int elmt=0;
  for (int i=0 ; i<positions.n ; i++) {

    zi = zj = positions.t[i][2];
    lay = ml.getMedium(zi);
    calcQintegralsSelf();
    for (int a=0 ; a<nbrQ ; a++) {
      Qintegrals.t[a][elmt] = Q1.t[a];
      Qintegrals.t[a][elmt+1] = Q2.t[a];
      Qintegrals.t[a][elmt+2] = Q3.t[a];
      Qintegrals.t[a][elmt+3] = Q4.t[a];
    }

    elmt+=4;
  }
  if (positions.n>1) {
    for (int i=0 ; i<positions.n ; i++) {
      for (int j=i+1 ; j<positions.n ; j++) {
        //cout << "Refl\t" << i << "\t" << j << endl;
        layI = ml.getMedium(positions.t[i][2]);
        layJ = ml.getMedium(positions.t[j][2]);
        dx = positions.t[i][0]-positions.t[j][0];
        dy = positions.t[i][1]-positions.t[j][1];
        zi = positions.t[i][2];
        zj = positions.t[j][2];
        bessj = Matrix(2*N+1,2*Ng);
        rauij = sqrt(dx*dx+dy*dy);
        for (int kk=0 ; kk<Ng ; kk++) {
          bes.assign(2*N,rauij*kappa.t[kk]);
          for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
        }
        for (int kk=Ng ; kk<2*Ng ; kk++) {
          bes.assign(2*N,rauij*kappa.t[kk]);
          for (int n=0 ; n<=2*N ; n++) bessj.t[n][kk] = bes.J.t[n];
        }
        calcQintegrals();
        for (int a=0 ; a<nbrQ ; a++) {
          Qintegrals.t[a][elmt] = Q1.t[a];
          Qintegrals.t[a][elmt+1] = Q2.t[a];
          Qintegrals.t[a][elmt+2] = Q3.t[a];
          Qintegrals.t[a][elmt+3] = Q4.t[a];
          //Qintegrals.t[a][elmt+1] = Qi.t[a];
        }
        elmt+=4;
      }
    }
  }
  return;
}

void Reflection1::calcRMatrices(const int &i, const int &j) {
  Q1 = Q2 = Q3 = Q4 = Vector(nbrQ,0.);
  if (i==j) {
    for (int a=0 ; a<nbrQ ; a++) {
      Q1.t[a] = Qintegrals.t[a][i*4];
      Q2.t[a] = Qintegrals.t[a][i*4+1];
      Q3.t[a] = Qintegrals.t[a][i*4+2];
      Q4.t[a] = Qintegrals.t[a][i*4+3];
    }
    calcRSelf();
  }
  else {
    dx = positions.t[i][0]-positions.t[j][0];
    dy = positions.t[i][1]-positions.t[j][1];
    rauij = sqrt(dx*dx+dy*dy);
    eiphi = (dx+i_*dy)/rauij;
    if (rauij==0.) eiphi = 1.;
    int elmt=getRefl(i,j);
    for (int a=0 ; a<nbrQ ; a++) {
      Q1.t[a] = Qintegrals.t[a][elmt];
      Q2.t[a] = Qintegrals.t[a][elmt+1];
      Q3.t[a] = Qintegrals.t[a][elmt+2];
      Q4.t[a] = Qintegrals.t[a][elmt+3];
    }
    calcR();
  }
}

Complex Reflection1::getQelmt(const int &u, const int &v, const int &w, const int &pt, const int &pm) {
  if (w<abs(u) || w<v) return 0.;
  int i=0,elmt;
  for (int U=-2 ; U<=2 ; U+=2) {
    for (int W=abs(U) ; W<=2*N ; W++) {
      for (int V=0 ; V<=W ; V++) {
        if(U==u && V==v && W==w) elmt = i;
        i++;
      }
    }
  }
  if ((pt%2==0) && (pm%2 == 0)) return Q1.t[elmt];
  else if ((pt%2==0) && (pm%2 != 0)) return Q2.t[elmt];
  else if ((pt%2!=0) && (pm%2 == 0)) return Q3.t[elmt];
  else if ((pt%2!=0) && (pm%2 != 0)) return Q4.t[elmt];
}

void Reflection1::calcQintegrals() {
  Q1 = Q2 = Q3 = Q4 = Vector(nbrQ,0.);
  Complex mult,mult1;
  Vector a1(2,1.),a2(2,1.);
  a2.t[1] = -1.;
  Vector ap1,as1,ap2,as2;
  int elmt;
  for (int k=0 ; k<2*Ng ; k++) {
    ml.setKxy(kappa.t[k]);
    ap1 = ml.sourceInt(1,a1,zj,zi);
    as1 = ml.sourceInt(0,a1,zj,zi);
    ap2 = ml.sourceInt(1,a2,zj,zi);
    as2 = ml.sourceInt(0,a2,zj,zi);

    elmt = 0;
    mult1 = gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[k+2*layJ*Ng];
    for (int u=-2 ; u<=2 ; u+=2) {
      for (int w=abs(u) ; w<=2*N ; w++) {
        for (int v=0 ; v<=w ; v++) {
          mult = mult1*powi(abs(v))*bessj.t[int(abs(v))][k]*duvw.t[k+2*layJ*Ng][elmt];
          Q1.t[elmt] += mult*((ap1.t[0]+powi(u)*as2.t[0]) + (-ap1.t[1]+powi(u)*as2.t[1]));
          Q2.t[elmt] += mult*((ap1.t[0]+powi(u)*as2.t[0]) - (-ap1.t[1]+powi(u)*as2.t[1]));
          Q3.t[elmt] += mult*((ap2.t[0]+powi(u)*as1.t[0]) + (-ap2.t[1]+powi(u)*as1.t[1]));
          Q4.t[elmt] += mult*((ap2.t[0]+powi(u)*as1.t[0]) - (-ap2.t[1]+powi(u)*as1.t[1]));
          elmt++;
        }
      }
    }
  }
  return;
}

void Reflection1::calcQintegralsSelf() {
  Q1 = Q2 = Q3 = Q4 = Vector(nbrQ,0.);
  Complex mult,mult1;
  Vector a1(2,1.),a2(2,1.);
  a2.t[1] = -1.;
  Vector ap1,as1,ap2,as2;
  int elmt;
  for (int k=0 ; k<2*Ng ; k++) {

    ml.setKxy(kappa.t[k]);
    ap1 = ml.sourceInt(1,a1,zj,zj);
    as1 = ml.sourceInt(0,a1,zj,zj);
    ap2 = ml.sourceInt(1,a2,zj,zj);
    as2 = ml.sourceInt(0,a2,zj,zj);
    elmt = 0;
    mult1 = gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[lay]/kz.t[k+2*lay*Ng];
    for (int u=-2 ; u<=2 ; u+=2) {
      for (int w=abs(u) ; w<=2*N ; w++) {
        for (int v=0 ; v<=w ; v++) {
          if (v==0) {
            mult = mult1*duvw.t[k+2*lay*Ng][elmt];
            Q1.t[elmt] += mult*((ap1.t[0]+powi(u)*as2.t[0]) + (-ap1.t[1]+powi(u)*as2.t[1]));
            Q2.t[elmt] += mult*((ap1.t[0]+powi(u)*as2.t[0]) - (-ap1.t[1]+powi(u)*as2.t[1]));
            Q3.t[elmt] += mult*((ap2.t[0]+powi(u)*as1.t[0]) + (-ap2.t[1]+powi(u)*as1.t[1]));
            Q4.t[elmt] += mult*((ap2.t[0]+powi(u)*as1.t[0]) - (-ap2.t[1]+powi(u)*as1.t[1]));
          }
          elmt++;
        }
      }
    }
  }
  return;
}

void Reflection1::calcRSelf() {
  R = Matrix(2*NM,2*NM,0.);
  int line=0,col=0;
  Complex mult,multi,Qm2,Q0,Qp2,Qm2i,Q0i,Qp2i;
  Complex Qam2,Qa0,Qap2,Qbm2,Qb0,Qbp2,Qcm2,Qc0,Qcp2,Qdm2,Qd0,Qdp2;
  double C1,C2,pm1;
  line = 0;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      col=0;
      for (int l=1 ; l<=N ; l++) {
        for (int k=-l ; k<=l ; k++) {
          mult = 0.25*powm1(m)*powi(l-n)*sqrt(double(n+n+1)*double(l+l+1));
          if (m==k) {
            for (int w=abs(n-l) ; w<=n+l ; w++) {
              C1 = cg.coef(w,n,m,l,-k)*cg.coef(w,n,1,l,-1);
              C2 = cg.coef(w,n,m,l,-k)*cg.coef(w,n,1,l,1);
              pm1 = powm1(n+l+w);
              Qam2 = getQelmt(-2,0,w,n-m,l-k+1);
              Qa0 = getQelmt(0,0,w,n-m,l-k+1);
              Qap2 = getQelmt(2,0,w,n-m,l-k+1);

              Qbm2 = getQelmt(-2,0,w,n-m,l-k);
              Qb0 = getQelmt(0,0,w,n-m,l-k);
              Qbp2 = getQelmt(2,0,w,n-m,l-k);

              Qcm2 = getQelmt(-2,0,w,n-m+1,l-k+1);
              Qc0 = getQelmt(0,0,w,n-m+1,l-k+1);
              Qcp2 = getQelmt(2,0,w,n-m+1,l-k+1);

              Qdm2 = getQelmt(-2,0,w,n-m+1,l-k);
              Qd0 = getQelmt(0,0,w,n-m+1,l-k);
              Qdp2 = getQelmt(2,0,w,n-m+1,l-k);

              R.t[line][col] += mult*(C1*(1.+pm1)*Qa0 + C2*(Qap2 + pm1*Qam2));
              R.t[line][col+NM] += mult*(C1*(1.-pm1)*Qb0 - C2*(Qbp2 - pm1*Qbm2));
              R.t[line+NM][col] += mult*(C1*(1.-pm1)*Qc0 + C2*(Qcp2 - pm1*Qcm2));
              R.t[line+NM][col+NM] += mult*(C1*(1.+pm1)*Qd0 - C2*(Qdp2 + pm1*Qdm2));
            }
          }

          col++;
        }
      }
      line++;
    }
  }
  return;
}

void Reflection1::calcR() {
  R = Ri = Matrix(2*NM,2*NM,0.);
  int line=0,col=0;
  Complex mult,multi,Qm2,Q0,Qp2,Qm2i,Q0i,Qp2i;
  Complex Qam2,Qa0,Qap2,Qbm2,Qb0,Qbp2,Qcm2,Qc0,Qcp2,Qdm2,Qd0,Qdp2;
  double C1,C2,pm1;
  line = 0;
  for (int n=1 ; n<=N ; n++) {
    for (int m=-n ; m<=n ; m++) {
      col=0;
      for (int l=1 ; l<=N ; l++) {
        for (int k=-l ; k<=l ; k++) {
          mult = 0.25*powm1(k+1)*powi(l-n)*pow(eiphi,double(m-k))*sqrt(double(n+n+1)*double(l+l+1));
          //multi = 0.25*powm1(m+1)*powi(l-n)*eiphii*sqrt(double(n+n+1)*double(l+l+1));
            for (int w=abs(n-l) ; w<=n+l ; w++) {
              C1 = cg.coef(w,n,m,l,-k)*cg.coef(w,n,1,l,-1);
              C2 = cg.coef(w,n,m,l,-k)*cg.coef(w,n,1,l,1);
              pm1 = powm1(n+l+w);
              if (m-k<0) {
                Qam2 = powm1(m-k)*getQelmt(2,k-m,w,n-m,l-k+1);
                Qa0 = powm1(m-k)*getQelmt(0,k-m,w,n-m,l-k+1);
                Qap2 = powm1(m-k)*getQelmt(-2,k-m,w,n-m,l-k+1);

                Qbm2 = powm1(m-k)*getQelmt(2,k-m,w,n-m,l-k);
                Qb0 = powm1(m-k)*getQelmt(0,k-m,w,n-m,l-k);
                Qbp2 = powm1(m-k)*getQelmt(-2,k-m,w,n-m,l-k);

                Qcm2 = powm1(m-k)*getQelmt(2,k-m,w,n-m+1,l-k+1);
                Qc0 = powm1(m-k)*getQelmt(0,k-m,w,n-m+1,l-k+1);
                Qcp2 = powm1(m-k)*getQelmt(-2,k-m,w,n-m+1,l-k+1);

                Qdm2 = powm1(m-k)*getQelmt(2,k-m,w,n-m+1,l-k);
                Qd0 = powm1(m-k)*getQelmt(0,k-m,w,n-m+1,l-k);
                Qdp2 = powm1(m-k)*getQelmt(-2,k-m,w,n-m+1,l-k);

              }
              else {
                Qam2 = getQelmt(-2,m-k,w,n-m,l-k+1);
                Qa0 = getQelmt(0,m-k,w,n-m,l-k+1);
                Qap2 = getQelmt(2,m-k,w,n-m,l-k+1);

                Qbm2 = getQelmt(-2,m-k,w,n-m,l-k);
                Qb0 = getQelmt(0,m-k,w,n-m,l-k);
                Qbp2 = getQelmt(2,m-k,w,n-m,l-k);

                Qcm2 = getQelmt(-2,m-k,w,n-m+1,l-k+1);
                Qc0 = getQelmt(0,m-k,w,n-m+1,l-k+1);
                Qcp2 = getQelmt(2,m-k,w,n-m+1,l-k+1);

                Qdm2 = getQelmt(-2,m-k,w,n-m+1,l-k);
                Qd0 = getQelmt(0,m-k,w,n-m+1,l-k);
                Qdp2 = getQelmt(2,m-k,w,n-m+1,l-k);
              }
              R.t[line][col] += mult*(C1*(1.+pm1)*Qa0 + C2*(Qap2 + pm1*Qam2));
              R.t[line][col+NM] += mult*(C1*(1.-pm1)*Qb0 - C2*(Qbp2 - pm1*Qbm2));
              R.t[line+NM][col] += mult*(C1*(1.-pm1)*Qc0 + C2*(Qcp2 - pm1*Qcm2));
              R.t[line+NM][col+NM] += mult*(C1*(1.+pm1)*Qd0 - C2*(Qdp2 + pm1*Qdm2));

              Ri.t[line][col] += powm1(m+k)*mult*(C1*(1.+pm1)*Qa0 + C2*(Qap2 + pm1*Qam2));
              Ri.t[line][col+NM] += powm1(m+k)*mult*(C1*(1.-pm1)*Qb0 - C2*(Qbp2 - pm1*Qbm2));
              Ri.t[line+NM][col] += powm1(m+k)*mult*(C1*(1.-pm1)*Qc0 + C2*(Qcp2 - pm1*Qcm2));
              Ri.t[line+NM][col+NM] += powm1(m+k)*mult*(C1*(1.+pm1)*Qd0 - C2*(Qdp2 + pm1*Qdm2));
          }
          col++;
        }
      }
      line++;
    }
  }
  return;
}

void Reflection1::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector app,asp,apt,ast,a(2,1.);
  Complex rpp,rpt,rsp,rst;
  Complex cost,mult,multExp,multExp1,multExp2;
  double alpha = atan2(dy,dx);
  int elmtJ,elmtI;
  Complex cosb;
  rauij = sqrt(dx*dx+dy*dy);
  eiphi = (dx+i_*dy)/rauij;
  if (rauij==0.) eiphi=1.;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  //ml.setLambda(lambda);
  for (int k=0 ; k<2*Ng ; k++) {

    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    bes.assign(N+1,rauij*kappa.t[k]);
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[elmtI]/ml.ki.t[layI];
    int i=0;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        a.t[0] = 1.;
        a.t[1] = powm1(n-m);
        app = ml.sourceInt(1,a,zzj,zzi);
        asp = ml.sourceInt(0,a,zzj,zzi);
        a.t[0] = 1.;
        a.t[1] = -powm1(n-m);
        apt = ml.sourceInt(1,a,zzj,zzi);
        ast = ml.sourceInt(0,a,zzj,zzi);

        multExp = mult*bes.J.t[int(abs(m))]*powi(-n+abs(m)-1)*pow(eiphi,double(m));
        multExp1 = mult*bes.J.t[int(abs(m+1))]*powi(-n+abs(m+1)-1)*pow(eiphi,double(m+1));
        multExp2 = mult*bes.J.t[int(abs(m-1))]*powi(-n+abs(m-1)-1)*pow(eiphi,double(m-1));

        Mmn.t[i][0] += 0.5*((multExp1+multExp2)*legPi.t[elmtJ][i]*(app.t[0]-app.t[1])*cosb - (multExp1-multExp2)*legTau.t[elmtJ][i]*(ast.t[0]+ast.t[1]));
        Mmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*legPi.t[elmtJ][i]*(app.t[0]-app.t[1])*cosb + (multExp1+multExp2)*i_*legTau.t[elmtJ][i]*(ast.t[0]+ast.t[1]));
        Mmn.t[i][2] += multExp*((app.t[0]+app.t[1])*legPi.t[elmtJ][i]*sqrt(1.-cosb*cosb));

        Nmn.t[i][0] += 0.5*((multExp1+multExp2)*legTau.t[elmtJ][i]*(apt.t[0]-apt.t[1])*cosb - (multExp1-multExp2)*legPi.t[elmtJ][i]*(asp.t[0]+asp.t[1]));
        Nmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*legTau.t[elmtJ][i]*(apt.t[0]-apt.t[1])*cosb + (multExp1+multExp2)*i_*legPi.t[elmtJ][i]*(asp.t[0]+asp.t[1]));
        Nmn.t[i][2] += multExp*(legTau.t[elmtJ][i]*(apt.t[0]+apt.t[1])*sqrt(1.-cosb*cosb));
        i++;
      }
    }
  }
  return;
}

/*void Reflection1::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector as,ap;
  Complex mult,multa,multExp,multExp1,multExp2;
  double alpha;
  int elmtJ,elmtI;
  Complex cosb;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  for (int k=0 ; k<2*Ng ; k++) {
    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    ap = ml.sourceInt(1,Vector(2,1.),zzj,zzi);
    as = ml.sourceInt(0,Vector(2,1.),zzj,zzi);
    //ap.t[0] = as.t[0] = exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
    //ap.t[1] = as.t[1] = 0.;
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[k+2*layI*Ng]/ml.ki.t[layI];
    for (int kk=0 ; kk<Ng ; kk++) {
      alpha = gaussLegAlpha.t[kk][0];
      multa = gaussLegAlpha.t[kk][1]*exp(i_*kappa.t[k]*(dx*cos(alpha)+dy*sin(alpha)));
      int i=0;
      for (int n=1 ; n<=N ; n++) {
        for (int m=-n ; m<=n ; m++) {
          multExp = mult*multa*powi(-n-1)*exp(i_*double(m)*alpha)*0.5/piDbl/sqrt(double(2*n*(n+1)));

          Mmn.t[i][0] += multExp*(double(m)*legmPi.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*cosb*cos(alpha) - i_*legTau.t[elmtJ][i]*(as.t[0]-powm1(n-m)*as.t[1])*sin(alpha));
          Mmn.t[i][1] += multExp*(double(m)*legmPi.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*cosb*sin(alpha) + i_*legTau.t[elmtJ][i]*(as.t[0]-powm1(n-m)*as.t[1])*cos(alpha));
          Mmn.t[i][2] -= multExp*(double(m)*legmPi.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*sqrt(1.-cosb*cosb));

          Nmn.t[i][0] += multExp*(legTau.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*cosb*cos(alpha) - i_*double(m)*legmPi.t[elmtJ][i]*(as.t[0]+powm1(n-m)*as.t[1])*sin(alpha));
          Nmn.t[i][1] += multExp*(legTau.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*cosb*sin(alpha) + i_*double(m)*legmPi.t[elmtJ][i]*(as.t[0]+powm1(n-m)*as.t[1])*cos(alpha));
          Nmn.t[i][2] -= multExp*(legTau.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*sqrt(1.-cosb*cosb));
          i++;
        }
      }
    }

  }
  return;
}*/

/*void Reflection1::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector as,ap;
  Complex cost,mult,multExp,multExp1,multExp2,app,apm,asp,asm;
  double alpha = atan2(dy,dx);
  int elmtJ,elmtI;
  Complex cosb;
  rauij = sqrt(dx*dx+dy*dy);
  eiphi = (dx+i_*dy)/rauij;
  if (rauij==0.) eiphi=1.;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  ml.setLambda(lambda);
  for (int k=0 ; k<2*Ng ; k++) {
    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    ap = ml.sourceInt(1,Vector(2,1.),zzj,zzi);
    as = ml.sourceInt(0,Vector(2,1.),zzj,zzi);

    bes.assign(N+1,rauij*kappa.t[k]);
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[elmtI]/ml.ki.t[layI];
    int i=0;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        app = powm1(n-m)*ap.t[0]+ap.t[1];
        apm = powm1(n-m)*ap.t[0]-ap.t[1];
        asp = powm1(n-m)*as.t[0]+as.t[1];
        asm = powm1(n-m)*as.t[0]-as.t[1];
        multExp = mult*bes.J.t[int(abs(m))]*powi(-n+abs(m)-1)*pow(eiphi,double(m))/sqrt(double(2*n*(n+1)));
        multExp1 = mult*bes.J.t[int(abs(m+1))]*powi(-n+abs(m+1)-1)*pow(eiphi,double(m+1))/sqrt(double(2*n*(n+1)));
        multExp2 = mult*bes.J.t[int(abs(m-1))]*powi(-n+abs(m-1)-1)*pow(eiphi,double(m-1))/sqrt(double(2*n*(n+1)));
        Mmn.t[i][0] += 0.5*((multExp1+multExp2)*double(m)*legmPi.t[elmtJ][i]*app*cosb + (multExp1-multExp2)*legTau.t[elmtJ][i]*asp);
        Mmn.t[i][1] -= 0.5*i_*((multExp1-multExp2)*double(m)*legmPi.t[elmtJ][i]*app*cosb + (multExp1+multExp2)*legTau.t[elmtJ][i]*asp);
        Mmn.t[i][2] += multExp*(double(m)*legmPi.t[elmtJ][i]*apm*sqrt(1.-cosb*cosb));

        Nmn.t[i][0] -= 0.5*((multExp1+multExp2)*legTau.t[elmtJ][i]*apm*cosb + (multExp1-multExp2)*double(m)*legmPi.t[elmtJ][i]*asm);
        Nmn.t[i][1] += 0.5*i_*((multExp1-multExp2)/i_*legTau.t[elmtJ][i]*apm*cosb - (multExp1+multExp2)*i_*double(m)*legmPi.t[elmtJ][i]*asm);
        Nmn.t[i][2] -= multExp*(legTau.t[elmtJ][i]*app*sqrt(1.-cosb*cosb));
        i++;
      }
    }
  }
  return;
}*/

/*void Reflection1::calcVSWFsIntegral(const double &dx, const double &dy, const double &zzi, const double &zzj) {
  Mmn = Nmn = Matrix(NM,3,0.);
  Vector as,ap;
  Complex cost,mult,multExp,multExp1,multExp2;
  double alpha = atan2(dy,dx);
  int elmtJ,elmtI;
  Complex cosb;
  rauij = sqrt(dx*dx+dy*dy);
  eiphi = (dx+i_*dy)/rauij;
  if (rauij==0.) eiphi=1.;
  layI = ml.getMedium(zzi);
  layJ = ml.getMedium(zzj);
  ml.setLambda(lambda);
  for (int k=0 ; k<2*Ng ; k++) {
    ml.setKxy(kappa.t[k]);
    elmtJ = k+2*layJ*Ng;
    elmtI = k+2*layI*Ng;
    ap = ml.sourceInt(1,Vector(2,1.),zzj,zzi);
    as = ml.sourceInt(0,Vector(2,1.),zzj,zzi);
    /*if (layI==layJ) {
      if (zzj<zzi) {
        ap.t[0] += exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
        as.t[0] += exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
      }
      else if (zzj>zzi) {
        ap.t[1] += exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
        as.t[1] += exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
      }
    }
    //ap.t[1] = as.t[1] = exp(i_*kz.t[k+2*layI*Ng]*fabs(zzi-zzj));
    //ap.t[0] = as.t[0] = 0.;
    bes.assign(N+1,rauij*kappa.t[k]);
    mult = -gaussLeg.t[k][1]*kappa.t[k]*dkappa.t[k]/ml.ki.t[layJ]/kz.t[elmtJ];
    cosb = kz.t[elmtI]/ml.ki.t[layI];
    int i=0;
    //ap.t[1] = as.t[1] = 0.;
    for (int n=1 ; n<=N ; n++) {
      for (int m=-n ; m<=n ; m++) {
        multExp = mult*bes.J.t[int(abs(m))]*powi(-n+abs(m)-1)*pow(eiphi,double(m))/sqrt(double(2*n*(n+1)));
        multExp1 = mult*bes.J.t[int(abs(m+1))]*powi(-n+abs(m+1)-1)*pow(eiphi,double(m+1))/sqrt(double(2*n*(n+1)));
        multExp2 = mult*bes.J.t[int(abs(m-1))]*powi(-n+abs(m-1)-1)*pow(eiphi,double(m-1))/sqrt(double(2*n*(n+1)));
        Mmn.t[i][0] += 0.5*((multExp1+multExp2)*double(m)*legmPi.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*cosb - (multExp1-multExp2)*legTau.t[elmtJ][i]*(as.t[0]-powm1(n-m)*as.t[1]));
        Mmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*double(m)*legmPi.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*cosb + (multExp1+multExp2)*i_*legTau.t[elmtJ][i]*(as.t[0]-powm1(n-m)*as.t[1]));
        Mmn.t[i][2] -= multExp*(double(m)*legmPi.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*sqrt(1.-cosb*cosb));

        Nmn.t[i][0] += 0.5*((multExp1+multExp2)*legTau.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*cosb - (multExp1-multExp2)*double(m)*legmPi.t[elmtJ][i]*(as.t[0]+powm1(n-m)*as.t[1]));
        Nmn.t[i][1] += 0.5*((multExp1-multExp2)/i_*legTau.t[elmtJ][i]*(ap.t[0]+powm1(n-m)*ap.t[1])*cosb + (multExp1+multExp2)*i_*double(m)*legmPi.t[elmtJ][i]*(as.t[0]+powm1(n-m)*as.t[1]));
        Nmn.t[i][2] -= multExp*(legTau.t[elmtJ][i]*(ap.t[0]-powm1(n-m)*ap.t[1])*sqrt(1.-cosb*cosb));
        i++;
      }
    }
  }
  return;
}*/
