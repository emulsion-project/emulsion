#include "Routines.h"

void calcCrossSections(const string &fname, const double &lambda1, const double &lambda2, const double &step,const bool &normalize, const double &cstNorm) {
  MultipleScattering ms(fname);
  ofstream of("out.txt");
  RVector Cexti,Cabsi,CextiDM,CextiQM,CextiSM,CextiDE,CextiQE,CextiSE;
  double Cabs,Cext,Csca;
  of << "lambda (nm) \t Cabs \t Csca \t Cext \n";
  for (double ll=lambda1 ; ll<=lambda2*1.0001 ; ll+=step) {
    ms.setLambda(ll);
    ms.solve();
    Cabs = Cext = 0.;
    Cexti = ms.calcCext(normalize,cstNorm);
    //Csca = ms.calcCsca(normalize,cstNorm);
    Cabsi = ms.calcCabs(normalize,cstNorm);
    for (int i=0 ; i<Cexti.n ; i++) Cext += Cexti.t[i];
    for (int i=0 ; i<Cabsi.n ; i++) Cabs += Cabsi.t[i];
    cout << ll*1e9 << "\t" << Cabs << "\t" << Cext-Cabs << "\t" << Cext << endl;
    of << 2.*piDbl*c*hbar/(e*ll) << "\t" << Cabs << "\t" << Cext-Cabs << "\t" << Cext << endl;
  }
}

void calcCrossSections(const string &fname, const double &lambda1, const double &lambda2, const double &step,const bool &normalize, const double &cstNorm, const int &nn, const int &pol) {
  MultipleScattering ms(fname);
  ofstream of("out.txt");
  RVector Cexti,Cabsi,CextiDM,CextiQM,CextiSM,CextiDE,CextiQE,CextiSE;
  double Cabs,Cext,Csca;
  of << "lambda (nm) \t Cabs \t Csca \t Cext \n";
  for (double ll=lambda1 ; ll<=lambda2*1.0001 ; ll+=step) {
    ms.setLambda(ll);
    ms.solve();
    Cabs = Cext = 0.;
    Cexti = ms.calcCextN(normalize,cstNorm,nn,pol);
    //Csca = ms.calcCsca(normalize,cstNorm);
    Cabsi = ms.calcCabsN(normalize,cstNorm,nn,pol);
    for (int i=0 ; i<Cexti.n ; i++) Cext += Cexti.t[i];
    for (int i=0 ; i<Cabsi.n ; i++) Cabs += Cabsi.t[i];
    cout << ll*1e9 << "\t" << Cabs << "\t" << Cext-Cabs << "\t" << Cext << endl;
    of << ll*1e9 << "\t" << Cabs << "\t" << Cext-Cabs << "\t" << Cext << endl;
  }
}

void calcComplexCext(const string &fname, const double &lambda1, const double &lambda2, const double &step,const bool &normalize, const double &cstNorm) {
  MultipleScattering ms(fname);
  ofstream of("out1.txt");
  Vector Cexti;
  of << "lambda (nm)";
  for (int i=0 ; i<ms.Np ; i++) {
    of << "\t" << "Cext" << i << "re" << "\t" << "Cext" << i << "im";
  }
  of << endl;
  for (double ll=lambda1 ; ll<=lambda2*1.0001 ; ll+=step) {
    ms.setLambda(ll);
    ms.solve();
    Cexti = ms.calcComplexCext(normalize,cstNorm);
    cout << ll*1e9;
    of << ll*1e9;
    for (int i=0 ; i<Cexti.n ; i++) {
      cout << "\t" << Cexti.t[i].re << "\t" << Cexti.t[i].im;
      of << "\t" << Cexti.t[i].re << "\t" << Cexti.t[i].im;
    }
    cout << endl;
    of << endl;
  }
}

void calcComplexCextN(const string &fname, const double &lambda1, const double &lambda2, const double &step,const bool &normalize, const double &cstNorm, const int &nn, const int &pol) {
  MultipleScattering ms(fname);
  ofstream of("out1.txt");
  Vector Cexti;
  of << "lambda (nm)";
  for (int i=0 ; i<ms.Np ; i++) {
    of << "\t" << "Cext" << i << "re" << "\t" << "Cext" << i << "im";
  }
  of << endl;
  Complex cexti;
  for (double ll=lambda1 ; ll<=lambda2*1.0001 ; ll+=step) {
    ms.setLambda(ll);
    ms.solve();
    Cexti = ms.calcComplexCextN(normalize,cstNorm,nn,pol);
    cexti = 0.;
    cout << ll*1e9;
    of << ll*1e9;
    for (int i=0 ; i<Cexti.n ; i++) {
      cexti += Cexti.t[i];
    }
    cout << "\t" << cexti.re << "\t" << cexti.im << endl;
    of << "\t" << cexti.re << "\t" << cexti.im << endl;
  }
}

void calcCabsi(const string &fname, const double &lambda1, const double &lambda2, const double &step,const bool &normalize, const double &cstNorm) {
  MultipleScattering ms(fname);
  ofstream of("out2.txt");
  RVector Cabsi;
  of << "lambda (nm)";
  for (int i=0 ; i<ms.Np ; i++) {
    of << "\t" << "Cabs" << i;
  }
  of << endl;
  for (double ll=lambda1 ; ll<=lambda2*1.0001 ; ll+=step) {
    ms.setLambda(ll);
    ms.solve();
    Cabsi = ms.calcCabs(normalize,cstNorm);
    cout << ll*1e9;
    of << ll*1e9;
    for (int i=0 ; i<Cabsi.n ; i++) {
      cout << "\t" << Cabsi.t[i];
      of << "\t" << Cabsi.t[i];
    }
    cout << endl;
    of << endl;
  }
}

void calcCextVSn(const string &fname, const int &N1, const int &N2, const int &pas, const bool &normalize, const double &cstNorm) {
  MultipleScattering ms(fname);
  ofstream of("out.txt");
  RVector Cexti;
  double Cext;
  of << "N \t Cext \n";
  for (int n=N1 ; n<=N2 ; n+=pas) {
    ms.setN(n);
    ms.solve();
    Cext = 0.;
    Cexti = ms.calcCext(normalize,cstNorm);
    for (int i=0 ; i<Cexti.n ; i++) Cext += Cexti.t[i];
    cout << n << "\t" << Cext << endl;
    of << n << "\t" << Cext << endl;
  }
}


void calcEnergyLossProbability(const string &fname, const double &e1, const double &e2, const double &step) {
  MultipleScattering ms(fname);
  ofstream of("out.txt");
  double proba,CL;
  of << "E (eV) \t Wavelength (nm) \t Proba \n";
  for (double ee=e1 ; ee<=e2*1.0001 ; ee+=step) {
    ms.setElectronEnergy(ee);
    ms.solve();
    proba = ms.calcEEL();
    CL = ms.calcCL();
    Complex epsilon = ms.particles[0].n2;
    cout << ee << "\t" << ms.lambda*1e9 << "\t" << proba << "\t" << CL << endl;
    of << ee << "\t" << ms.lambda*1e9 << "\t" << (proba) << "\t" << CL << endl;
  }
}

void calcNearFieldPlane(const string &fname, const int &N1, const int &N2,const double &t1,const double &t2,const double &p1,const double &p2,const int &plane, const bool &cart) {
  cout << "Solving the linear system...";
  MultipleScattering ms(fname);

  //ms.setLambda(480e-9);
  //ms.setElectronEnergy(2.91);
  ms.solve();

  //ms.calcPmn();
  cout << "done" << endl;
  //if(ms.incidence.type == 1) ms.incidence.setElectron(ms.N,ms.incidence.E0,ms.incidence.E,0.,0.);
  RMatrix X(N1,N2),Y(N1,N2),Z(N1,N2);
  RMatrix E1r(N1,N2),E1i(N1,N2),E2r(N1,N2),E2i(N1,N2),E3r(N1,N2),E3i(N1,N2);
  ofstream of1("near_field/nbrPlots.txt");
  of1 << 1 << "\t" << N1 << "\t" << N2;
  double pas1 = fabs(t1-t2)/double(N1-1);
  double pas2 = fabs(p1-p2)/double(N2-1);
  double R,T,P,xx,yy,zz;
  Vector3d field;
  cout << "Compute the near field : " << endl;
  for (int t=0 ; t<N1 ; t++) {
    cout << N1-t << "  lines remaining" << endl;
    T = t1+double(t)*pas1;
    for (int p=0 ; p<N2 ; p++) {
      //cout << N2-p << "  columns remaining" << endl;
      P = p1+double(p)*pas2;
      if (plane==0) { xx = T; yy = P ; zz = 0.; }
      else if (plane==1) { xx = T; yy = 0. ; zz = P; }
      else { xx = 0.; yy = T ; zz = P; }
      field = ms.calcEsca(xx,yy,zz) + ms.calcEinc(xx,yy,zz);
      //cout << "\t" << dotProduct(field,field) << endl;

      if (!cart) field = toSphericalCoordinates(field,T,P);
      X.t[t][p] = xx;
      Y.t[t][p] = yy;
      Z.t[t][p] = zz;
      E1r.t[t][p] = real(field.t[0]);
      E1i.t[t][p] = imag(field.t[0]);
      E2r.t[t][p] = real(field.t[1]);
      E2i.t[t][p] = imag(field.t[1]);
      E3r.t[t][p] = real(field.t[2]);
      E3i.t[t][p] = imag(field.t[2]);
    }
  }
  //cout << (dotProduct(ms.calcEinc(60e-9,0.,0.),ms.calcEinc(60e-9,0.,0.))) << endl;
  ofstream if1("near_field/x.txt");
  if1 << X;
  ofstream if2("near_field/y.txt");
  if2 << Y;
  ofstream if3("near_field/z.txt");
  if3 << Z;

  ofstream if4("near_field/E1r.txt");
  if4 << E1r;
  ofstream if5("near_field/E1i.txt");
  if5 << E1i;
  ofstream if6("near_field/E2r.txt");
  if6 << E2r;
  ofstream if7("near_field/E2i.txt");
  if7 << E2i;
  ofstream if8("near_field/E3r.txt");
  if8 << E3r;
  ofstream if9("near_field/E3i.txt");
  if9 << E3i;
}

void calcNearFieldParticleSurface(const string &fname, const int &N1, const int &N2, const bool &cart) {
  MultipleScattering ms(fname);
  cout << "Solving the linear system...";
  ms.setLambda(720e-9);
  ms.solve();
  cout << "done" << endl;
  if(ms.incidence.type == 1) ms.incidence.setElectron(ms.N,ms.incidence.E0,ms.incidence.E,0.,0.);
  RMatrix X(N1,N2*ms.Np),Y(N1,N2*ms.Np),Z(N1,N2*ms.Np);
  RMatrix E1r(N1,N2*ms.Np),E1i(N1,N2*ms.Np),E2r(N1,N2*ms.Np),E2i(N1,N2*ms.Np),E3r(N1,N2*ms.Np),E3i(N1,N2*ms.Np);
  ofstream of1("near_field/nbrPlots.txt");
  of1 << ms.Np << "\t" << N1 << "\t" << N2;
  double pas1 = piDbl/double(N1-1);
  double pas2 = 2.*piDbl/double(N2-1);
  double R,T,P,xx,yy,zz;
  Vector3d field;
  cout << "Compute the near field : " << endl;
  for (int i=0 ; i<ms.Np ; i++) {
    cout << ms.Np-i << "  particles remaining" << endl;
      for (int t=0 ; t<N1 ; t++) {
      T = double(t)*pas1;
      R = ms.particles[i].calcR(cos(T));
      //cout << T << "\n" << R << endl;
      zz = R*cos(T)+ms.z[i];
      for (int p=0 ; p<N2 ; p++) {
        P = double(p)*pas2;
        xx = R*cos(P)*sin(T) + ms.x[i];
        yy = R*sin(P)*sin(T) + ms.y[i];
        field = ms.calcEsca(xx,yy,zz);
        if (!cart) field = toSphericalCoordinates(field,T,P);
        X.t[t][p+N2*i] = xx;
        Y.t[t][p+N2*i] = yy;
        Z.t[t][p+N2*i] = zz;
        E1r.t[t][p+N2*i] = real(field.t[0]);
        E1i.t[t][p+N2*i] = imag(field.t[0]);
        E2r.t[t][p+N2*i] = real(field.t[1]);
        E2i.t[t][p+N2*i] = imag(field.t[1]);
        E3r.t[t][p+N2*i] = real(field.t[2]);
        E3i.t[t][p+N2*i] = imag(field.t[2]);
      }
    }
  }


  ofstream if1("near_field/x.txt");
  if1 << X;
  ofstream if2("near_field/y.txt");
  if2 << Y;
  ofstream if3("near_field/z.txt");
  if3 << Z;

  ofstream if4("near_field/E1r.txt");
  if4 << E1r;
  ofstream if5("near_field/E1i.txt");
  if5 << E1i;
  ofstream if6("near_field/E2r.txt");
  if6 << E2r;
  ofstream if7("near_field/E2i.txt");
  if7 << E2i;
  ofstream if8("near_field/E3r.txt");
  if8 << E3r;
  ofstream if9("near_field/E3i.txt");
  if9 << E3i;
}

void calcFarFieldSurface(const string &fname, const int &N1, const int &N2) {
  MultipleScattering ms(fname);
  cout << "Solving the linear system...";
  ms.setLambda(488e-9);
  //ms.setElectronEnergy(2.78);
  ms.solve();

  //ms.calcPmn();
  cout << "done" << endl;
  if(ms.incidence.type == 1) ms.incidence.setElectron(ms.N,ms.incidence.E0,ms.incidence.E,0.,0.);
  RMatrix X(N1,N2),Y(N1,N2),Z(N1,N2);
  RMatrix E1r(N1,N2),E1i(N1,N2),E2r(N1,N2),E2i(N1,N2),E3r(N1,N2),E3i(N1,N2);
  ofstream of1("near_field/nbrPlots.txt");
  of1 << 1 << "\t" << N1 << "\t" << N2;
  double pas1 = piDbl/double(N1-1);
  double pas2 = 2.*piDbl/double(N2-1);
  double R,T,P,xx,yy,zz;
  Vector3d field;
  cout << "Compute the far field : " << endl;
  for (int t=0 ; t<N1 ; t++) {
      T = double(t)*pas1;
      R = 1.;
      zz = R*cos(T);
      cout << t << endl;
      for (int p=0 ; p<N2 ; p++) {
        P = double(p)*pas2;
        xx = R*cos(P)*sin(T);
        yy = R*sin(P)*sin(T);
        field = ms.calcEscaInf(T,P);
        //field = toCartesianCoordinates(field,T,P);
        X.t[t][p] = xx;
        Y.t[t][p] = yy;
        Z.t[t][p] = zz;
        E1r.t[t][p] = real(field.t[0]);
        //cout << t << "\t" << p << "\t" << field << endl;
        E1i.t[t][p] = imag(field.t[0]);
        E2r.t[t][p] = real(field.t[1]);
        E2i.t[t][p] = imag(field.t[1]);
        E3r.t[t][p] = real(field.t[2]);
        E3i.t[t][p] = imag(field.t[2]);
      }
    }


  ofstream if1("near_field/x.txt");
  if1 << X;
  ofstream if2("near_field/y.txt");
  if2 << Y;
  ofstream if3("near_field/z.txt");
  if3 << Z;

  ofstream if4("near_field/E1r.txt");
  if4 << E1r;
  ofstream if5("near_field/E1i.txt");
  if5 << E1i;
  ofstream if6("near_field/E2r.txt");
  if6 << E2r;
  ofstream if7("near_field/E2i.txt");
  if7 << E2i;
  ofstream if8("near_field/E3r.txt");
  if8 << E3r;
  ofstream if9("near_field/E3i.txt");
  if9 << E3i;
}

void calcSurfaceCharges(const string &fname, const int &N1, const int &N2) {
  MultipleScattering ms(fname);
  //ms.setLambda(900e-9);
  //ms.setElectronEnergy(2.36);
  cout << "Solving the linear system...";
  ms.solve();

  cout << "done" << endl;
  if(ms.incidence.type == 1) ms.incidence.setElectron(ms.N,ms.incidence.E0,ms.incidence.E,0.,0.);
  RMatrix X(N1,N2*ms.Np),Y(N1,N2*ms.Np),Z(N1,N2*ms.Np);
  RMatrix E1(N1,N2*ms.Np),E2(N1,N2*ms.Np);
  ofstream of1("near_field/nbrPlots.txt");
  of1 << ms.Np << "\t" << N1 << "\t" << N2;
  double pas1 = piDbl/double(N1-1);
  double pas2 = 2.*piDbl/double(N2-1);
  double R,T,P,xx,yy,zz;
  Vector3d fieldExt,fieldInt;
  cout << "Compute the near field : " << endl;
  for (int i=0 ; i<ms.Np ; i++) {
    cout << ms.Np-i << "  particles remaining" << endl;
      for (int t=0 ; t<N1 ; t++) {
      T = double(t)*pas1;
      R = ms.particles[i].calcR(cos(T));


      zz = R*cos(T)+ms.z[i];
      for (int p=0 ; p<N2 ; p++) {
        P = double(p)*pas2;
        xx = R*cos(P)*sin(T) + ms.x[i];
        yy = R*sin(P)*sin(T) + ms.y[i];
        fieldExt = ms.calcEsca(xx,yy,zz)+ms.calcEinc(xx,yy,zz);
        fieldInt = ms.calcEint(i,xx,yy,zz);
        fieldExt = toSphericalCoordinates(fieldExt,T,P);
        fieldInt = toSphericalCoordinates(fieldInt,T,P);
        X.t[t][p+N2*i] = xx;
        Y.t[t][p+N2*i] = yy;
        Z.t[t][p+N2*i] = zz;
        //E1.t[t][p+N2*i] = real(fieldExt.t[0]*eps0*(ms.particles[i].n1*ms.particles[i].n1)-fieldInt.t[0]*eps0*(ms.particles[i].n2*ms.particles[i].n2));
        E1.t[t][p+N2*i] = real(fieldExt.t[0]-fieldInt.t[0]);
      }
    }
  }


  ofstream if1("near_field/x.txt");
  if1 << X;
  ofstream if2("near_field/y.txt");
  if2 << Y;
  ofstream if3("near_field/z.txt");
  if3 << Z;

  ofstream if4("near_field/E1r.txt");
  if4 << E1;
}

void calcEELSMapping(const string &fname, const int &N1, const int &N2,const double &t1,const double &t2,const double &p1,const double &p2) {
  MultipleScattering ms(fname);
  double XX,YY;
  RMatrix X(N1,N2),Y(N1,N2),Z(N1,N2);
  RMatrix EEL(N1,N2);
  ofstream of1("eelMap/nbrPlots.txt");
  of1 << 1 << "\t" << N1 << "\t" << N2;
  double pas1 = fabs(t1-t2)/double(N1-1);
  double pas2 = fabs(p1-p2)/double(N2-1);
  double R,T,P,xx,yy;
  double eel;
  bool calc = true;
  for (int t=0 ; t<N1 ; t++) {
    cout << t << endl;
    T = t1+double(t)*pas1;
    for (int p=0 ; p<N2 ; p++) {
      cout << p << " ";
      calc = true;
      P = p1+double(p)*pas2;
      xx = T; yy = P;

      for (int i=0 ; i<ms.Np ; i++) {
        XX = xx-ms.x[i];
        YY = yy-ms.y[i];
        R = sqrt(XX*XX+YY*YY);
        if (R < (ms.particles[i].geom.t[0])*0.9999999) calc = false;
      }
      if (calc) {


        double ee = 2.;
        ms.setLambda(2.*piDbl*c*hbar/(ee*e));
        ms.incidence.x0 = xx;
        ms.incidence.y0 = yy;

        //ms.incidence.setElectron(ms.N,ms.incidence.E0,ee,xx,yy);
        //ms.setElectronEnergy(2.76);
        //cout << ms.incidence.x0 << "\t" << ms.incidence.y0 << endl;
        ms.solve();
        eel = ms.calcEEL();
      }
      else eel = 0.;
      X.t[t][p] = xx;
      Y.t[t][p] = yy;
      Z.t[t][p] = 0.;
      EEL.t[t][p] = eel;
    }
    cout << "\n";
  }
  ofstream if1("eelMap/x.txt");
  if1 << X;
  ofstream if2("eelMap/y.txt");
  if2 << Y;
  ofstream if3("eelMap/z.txt");
  if3 << Z;
  ofstream if4("eelMap/E1r.txt");
  if4 << EEL;
}

Matrix polesAmpExtraction(const string &fname, const int &nPoles, const int &nPoints, const double &l1, const double &l2, const int &part) {
  MultipleScattering ms(fname);
  double pas = (l2-l1)/double(nPoints-1);
  Vector lambdas(nPoints);
  Vector omegas(nPoints);
  Vector cext(nPoints);
  Vector CextC;
  Matrix PA;
  double ll;
  for (int i=0 ; i<nPoints ; i++) {
    double ll = l2-double(i)*pas;
    ll = l1+double(i)*pas;
    lambdas.t[i] = ll;
    omegas.t[i] = 2.*piDbl*c/ll;
    ms.setLambda(ll);
    ms.solve();
    CextC = ms.calcComplexCext(false,1e-9);
  }
  PA = polesAmp(cext,omegas,nPoles);
  return PA;
}

Matrix polesAmpExtraction(MultipleScattering &ms, const int &nPoles, const int &nPoints, const double &l1, const double &l2, const int &part) {
  double pas = (l2-l1)/double(nPoints-1);
  Vector lambdas(nPoints);
  Vector omegas(nPoints);
  Vector cext(nPoints);
  Vector CextC;
  Matrix PA;
  double ll;
  for (int i=0 ; i<nPoints ; i++) {
    double ll = l2-double(i)*pas;
    ll = l1+double(i)*pas;
    lambdas.t[i] = ll;
    omegas.t[i] = 2.*piDbl*c/ll;
    ms.setLambda(ll);
    ms.solve();
    CextC = ms.calcComplexCext(false,1e-9);
    cext.t[i] = 0.;
    //for (int ii=0 ; ii<ms.Np ; ii++) cext.t[i] += CextC.t[ii];

    if (part == 0) cext.t[i] = CextC.t[0];
    else cext.t[i] = CextC.t[2];
    //if (part == 0) cext.t[i] = CextC.t[0]+CextC.t[1]+CextC.t[2];
    //else cext.t[i] = CextC.t[3]+CextC.t[4]+CextC.t[5];
  }
  PA = polesAmp(cext,omegas,nPoles);
  //for (int i=0 ; i<nPoles ; i++) PA.t[i][0] = 2.*piDbl*c/PA.t[i][0];
  return PA;
}


Vector poleAmpExtraction(const string &fname, const int &nPoints, const double &l1, const double &l2, const int &part) {
  MultipleScattering ms(fname);
  double pas = (l2-l1)/double(nPoints-1);
  Vector lambdas(nPoints);
  Vector omegas(nPoints);
  Vector cext(nPoints);
  Vector CextC;
  Vector PA;
  double ll;
  for (int i=0 ; i<nPoints ; i++) {
    double ll = l1+double(i)*pas;
    lambdas.t[i] = ll;
    omegas.t[i] = 2.*piDbl*c/ll;
    ms.setLambda(ll);
    ms.solve();
    CextC = ms.calcComplexCext(false,1e-9);
    cext.t[i] = CextC.t[part];
  }
  PA = poleAmplitude(omegas,cext);
  //PA.t[0] = 2.*piDbl*c/PA.t[0];
  return PA;
}

Vector poleAmpExtraction(const string &fname, const double &z, const int &nPoints, const double &l1, const double &l2, const int &n, const int &m, const int &elec) {
  MultipleScattering ms(fname);
  //ms.z[0] = z;
  //ms.x[1] = -z;
  double pas = (l2-l1)/double(nPoints-1);
  Vector lambdas(nPoints);
  Vector omegas(nPoints);
  Vector cext(nPoints);
  Vector CextC;
  Vector PA,Pmn;
  double ll;
  for (int i=0 ; i<nPoints ; i++) {
    double ll = l1+double(i)*pas;
    lambdas.t[i] = ll;
    omegas.t[i] = 2.*piDbl*c/ll;
    ms.setLambda(ll);
    ms.solve();
    //ms.calcPmn();
    //CextC = ms.calcComplexCext(false,1e-9);
    //cext.t[i] = CextC.t[part];
    Pmn = ms.Pmni;
    ms.trans.setParameters(1,false,ms.N,ms.particles[0].k1,0.,0.,z);
    Pmn = ms.trans.translate(Pmn,0);
    cext.t[i] = Pmn.t[ii(n,m)+elec*ms.NM];
  }
  PA = poleAmplitude(omegas,cext);
  PA.t[0] = 2.*piDbl*c/PA.t[0];
  return PA;
}
