#include <iostream>
#include <fstream>
#include "Multilayer.h"

using namespace std;

int main()
{
  string str;
  ifstream ifs("multilayerConfig.txt");
  double lambda,angle,bottomPos;
  int nbrMediums;
  RVector Zi,thicknesses;
  IVector mediums;

  ifs >> str >> lambda >> str >> angle >> str >> nbrMediums >> str;
  angle *= piDbl/180.;
  int i;
  if (nbrMediums>2) {
    thicknesses = RVector(nbrMediums-2);
    for (i=0 ; i<nbrMediums-2 ; i++) ifs >> thicknesses.t[i];
  }
  else {
    thicknesses = RVector(1);
    ifs >> thicknesses.t[0];
  }

  if (nbrMediums>1) {
    Zi = RVector(nbrMediums-1);
    Zi.t[0] = 0.;
    if (nbrMediums>2) for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
  }
  else {
    Zi = RVector(1);
    Zi.t[0]=0.;
  }
  ifs >> str;
  mediums = IVector(nbrMediums);
  for (i=0 ; i<nbrMediums ; i++) ifs >> mediums.t[i];

  Multilayer ml(nbrMediums,Zi,mediums,lambda);

  int nbrVariables,var1,var2,var1layerthick;
  double var1min,var1max,var1step,var2min,var2max,var2step;
  ifs >> str >> nbrVariables;
  Vector ampl(2),aOut;
  int counting = 0;
  double totCalc;

  if (nbrVariables == 1) {
    ifs >> str >> var1;
    ifs >> var1min >> var1max >> var1step;
    if (var1 == 2) ifs >> var1layerthick;
    if (var1 == 0) {
      ofstream ofS("outputs/variable_wavelength/variable_wavelength_polS.txt");
      ofstream ofP("outputs/variable_wavelength/variable_wavelength_polP.txt");
      ofS << "Wavelength (nm) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      ofP << "Wavelength (nm) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      for (double ll=var1min ; ll<=var1max+var1step/100. ; ll+=var1step) {
        ml.setLambda(ll);
        if (angle<piDbl) {
          ml.setKxy(ml.ki.t[0]*sin(angle));
          ampl.t[0] = 1.; ampl.t[1] = 0.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll*1e9 << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll*1e9 << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
        }
        else {
          ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
          ampl.t[0] = 0.; ampl.t[1] = 1.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll*1e9 << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll*1e9 << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
        }
      }
    }
    else if (var1 == 1) {
      ofstream ofS("outputs/variable_angle/variable_angle_polS.txt");
      ofstream ofP("outputs/variable_angle/variable_angle_polP.txt");
      ofS << "Angle (°) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      ofP << "Angle (°) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      for (double ll=var1min ; ll<=var1max+var1step/100. ; ll+=var1step) {
        angle = ll*piDbl/180.;
        if (angle<piDblby2) {
          ml.setKxy(ml.ki.t[0]*sin(angle));
          ampl.t[0] = 1.; ampl.t[1] = 0.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
        }
        else {
          ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
          ampl.t[0] = 0.; ampl.t[1] = 1.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
        }
      }
    }
    else if (var1 == 2) {
      ofstream ofS("outputs/variable_thickness/variable_thickness_polS.txt");
      ofstream ofP("outputs/variable_thickness/variable_thickness_polP.txt");
      ofS << "Thickness (nm) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      ofP << "Thickness (nm) \t Re(r) \t Im(r) \t Re(t) \t Im(t) \t R \t T" << endl;
      for (double ll=var1min ; ll<=var1max+var1step/100. ; ll+=var1step) {
        thicknesses.t[var1layerthick-1] = ll;
        for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
        ml.setParameters(nbrMediums,Zi,mediums,lambda);
        if (angle<piDblby2) {
          ml.setKxy(ml.ki.t[0]*sin(angle));
          ampl.t[0] = 1.; ampl.t[1] = 0.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << norm(aOut.t[0]) << "\t" << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << endl;
        }
        else {
          ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
          ampl.t[0] = 0.; ampl.t[1] = 1.;
          aOut = ml.sourceExtOut(ampl,0);
          ofS << ll << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
          aOut = ml.sourceExtOut(ampl,1);
          ofP << ll << "\t" << aOut.t[1].re << "\t" << aOut.t[1].im << "\t" << aOut.t[0].re << "\t" << aOut.t[0].im << "\t" << norm(aOut.t[1]) << "\t" << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << endl;
        }
      }
    }
  }

  else if (nbrVariables == 2) {
    ifs >> str >> var1 >> var2;
    ifs >> var1min >> var1max >> var1step;
    if (var1 == 2) ifs >> var1layerthick;
    ifs >> var2min >> var2max >> var2step;
    if (var2 == 2) ifs >> var1layerthick;
    totCalc = ((var1max-var1min)/var1step+1.)*((var2max-var2min)/var2step+1.);
    counting = 0;

    if ((var1 == 0 && var2 == 1) || (var1 == 1 && var2 == 0)) {
      ofstream ofSrer("outputs/variable_wavelength_angle/re(r)_polS.txt");
      ofstream ofPrer("outputs/variable_wavelength_angle/re(r)_polP.txt");
      ofstream ofSimr("outputs/variable_wavelength_angle/im(r)_polS.txt");
      ofstream ofPimr("outputs/variable_wavelength_angle/im(r)_polP.txt");
      ofstream ofSret("outputs/variable_wavelength_angle/re(t)_polS.txt");
      ofstream ofPret("outputs/variable_wavelength_angle/re(t)_polP.txt");
      ofstream ofSimt("outputs/variable_wavelength_angle/im(t)_polS.txt");
      ofstream ofPimt("outputs/variable_wavelength_angle/im(t)_polP.txt");
      ofstream ofSR("outputs/variable_wavelength_angle/R_polS.txt");
      ofstream ofPR("outputs/variable_wavelength_angle/R_polP.txt");
      ofstream ofST("outputs/variable_wavelength_angle/T_polS.txt");
      ofstream ofPT("outputs/variable_wavelength_angle/T_polP.txt");
      ofstream ofWavelength("outputs/variable_wavelength_angle/wavelength.txt");
      ofstream ofAngle("outputs/variable_wavelength_angle/angle.txt");

      for (double ll1 = var1min ; ll1 <= var1max+var1step*0.01 ; ll1+=var1step) {
        for (double ll2 = var2min ; ll2 <= var2max+var2step*0.01 ; ll2+=var2step) {
          if (var1 == 0) {
            ml.setLambda(ll1);
            angle = ll2*piDbl/180.;
            ofWavelength << ll1*1e9 << "\t";
            ofAngle << ll2 << "\t";
          }
          else {
            ml.setLambda(ll2);
            angle = ll1*piDbl/180.;
            ofWavelength << ll2*1e9 << "\t";
            ofAngle << ll1 << "\t";
          }
          if (angle<piDblby2) {
            ml.setKxy(ml.ki.t[0]*sin(angle));
            ampl.t[0] = 1.; ampl.t[1] = 0.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[0].re << "\t";
            ofSimr << aOut.t[0].im << "\t";
            ofSret << aOut.t[1].re << "\t";
            ofSimt << aOut.t[1].im << "\t";
            ofSR << norm(aOut.t[0]) << "\t";
            ofST << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[0].re << "\t";
            ofPimr << aOut.t[0].im << "\t";
            ofPret << aOut.t[1].re << "\t";
            ofPimt << aOut.t[1].im << "\t";
            ofPR << norm(aOut.t[0]) << "\t";
            ofPT << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
          }
          else {
            ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
            ampl.t[0] = 0.; ampl.t[1] = 1.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[1].re << "\t";
            ofSimr << aOut.t[1].im << "\t";
            ofSret << aOut.t[0].re << "\t";
            ofSimt << aOut.t[0].im << "\t";
            ofSR << norm(aOut.t[1]) << "\t";
            ofST << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[1].re << "\t";
            ofPimr << aOut.t[1].im << "\t";
            ofPret << aOut.t[0].re << "\t";
            ofPimt << aOut.t[0].im << "\t";
            ofPR << norm(aOut.t[1]) << "\t";
            ofPT << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
          }
          counting++;
        }
        cout << double(counting)/totCalc*100. << " %" << endl;
        ofWavelength << endl;
        ofAngle << endl;
        ofSrer << endl;
        ofSimr << endl;
        ofSret << endl;
        ofSimt << endl;
        ofSR << endl;
        ofST << endl;
        ofPrer << endl;
        ofPimr << endl;
        ofPret << endl;
        ofPimt << endl;
        ofPR << endl;
        ofPT << endl;
      }
    }
    else if ((var1 == 0 && var2 == 2) || (var1 == 2 && var2 == 0)) {
      ofstream ofSrer("outputs/variable_wavelength_thickness/re(r)_polS.txt");
      ofstream ofPrer("outputs/variable_wavelength_thickness/re(r)_polP.txt");
      ofstream ofSimr("outputs/variable_wavelength_thickness/im(r)_polS.txt");
      ofstream ofPimr("outputs/variable_wavelength_thickness/im(r)_polP.txt");
      ofstream ofSret("outputs/variable_wavelength_thickness/re(t)_polS.txt");
      ofstream ofPret("outputs/variable_wavelength_thickness/re(t)_polP.txt");
      ofstream ofSimt("outputs/variable_wavelength_thickness/im(t)_polS.txt");
      ofstream ofPimt("outputs/variable_wavelength_thickness/im(t)_polP.txt");
      ofstream ofSR("outputs/variable_wavelength_thickness/R_polS.txt");
      ofstream ofPR("outputs/variable_wavelength_thickness/R_polP.txt");
      ofstream ofST("outputs/variable_wavelength_thickness/T_polS.txt");
      ofstream ofPT("outputs/variable_wavelength_thickness/T_polP.txt");
      ofstream ofWavelength("outputs/variable_wavelength_thickness/wavelength.txt");
      ofstream ofThickness("outputs/variable_wavelength_thickness/thickness.txt");

      for (double ll1 = var1min ; ll1 <= var1max+var1step*0.01 ; ll1+=var1step) {
        for (double ll2 = var2min ; ll2 <= var2max+var2step*0.01 ; ll2+=var2step) {
          if (var1 == 0) {
            ml.setLambda(ll1);
            thicknesses.t[var1layerthick-1] = ll2;
            for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
            ml.setParameters(nbrMediums,Zi,mediums,lambda);
            ofWavelength << ll1*1e9 << "\t";
            ofThickness << ll2*1e9 << "\t";
          }
          else {
            ml.setLambda(ll2);
            thicknesses.t[var1layerthick-1] = ll1;
            for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
            ml.setParameters(nbrMediums,Zi,mediums,lambda);
            ofWavelength << ll2*1e9 << "\t";
            ofThickness << ll1*1e9 << "\t";
          }
          if (angle<piDblby2) {
            ml.setKxy(ml.ki.t[0]*sin(angle));
            ampl.t[0] = 1.; ampl.t[1] = 0.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[0].re << "\t";
            ofSimr << aOut.t[0].im << "\t";
            ofSret << aOut.t[1].re << "\t";
            ofSimt << aOut.t[1].im << "\t";
            ofSR << norm(aOut.t[0]) << "\t";
            ofST << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[0].re << "\t";
            ofPimr << aOut.t[0].im << "\t";
            ofPret << aOut.t[1].re << "\t";
            ofPimt << aOut.t[1].im << "\t";
            ofPR << norm(aOut.t[0]) << "\t";
            ofPT << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
          }
          else {
            ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
            ampl.t[0] = 0.; ampl.t[1] = 1.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[1].re << "\t";
            ofSimr << aOut.t[1].im << "\t";
            ofSret << aOut.t[0].re << "\t";
            ofSimt << aOut.t[0].im << "\t";
            ofSR << norm(aOut.t[1]) << "\t";
            ofST << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[1].re << "\t";
            ofPimr << aOut.t[1].im << "\t";
            ofPret << aOut.t[0].re << "\t";
            ofPimt << aOut.t[0].im << "\t";
            ofPR << norm(aOut.t[1]) << "\t";
            ofPT << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
          }
          counting++;
        }
        cout << double(counting)/totCalc*100. << " %" << endl;
        ofWavelength << endl;
        ofThickness << endl;
        ofSrer << endl;
        ofSimr << endl;
        ofSret << endl;
        ofSimt << endl;
        ofSR << endl;
        ofST << endl;
        ofPrer << endl;
        ofPimr << endl;
        ofPret << endl;
        ofPimt << endl;
        ofPR << endl;
        ofPT << endl;
      }
    }
    else if ((var1 == 1 && var2 == 2) || (var1 == 2 && var2 == 1)) {
      ofstream ofSrer("outputs/variable_angle_thickness/re(r)_polS.txt");
      ofstream ofPrer("outputs/variable_angle_thickness/re(r)_polP.txt");
      ofstream ofSimr("outputs/variable_angle_thickness/im(r)_polS.txt");
      ofstream ofPimr("outputs/variable_angle_thickness/im(r)_polP.txt");
      ofstream ofSret("outputs/variable_angle_thickness/re(t)_polS.txt");
      ofstream ofPret("outputs/variable_angle_thickness/re(t)_polP.txt");
      ofstream ofSimt("outputs/variable_angle_thickness/im(t)_polS.txt");
      ofstream ofPimt("outputs/variable_angle_thickness/im(t)_polP.txt");
      ofstream ofSR("outputs/variable_angle_thickness/R_polS.txt");
      ofstream ofPR("outputs/variable_angle_thickness/R_polP.txt");
      ofstream ofST("outputs/variable_angle_thickness/T_polS.txt");
      ofstream ofPT("outputs/variable_angle_thickness/T_polP.txt");
      ofstream ofThickness("outputs/variable_angle_thickness/thickness.txt");
      ofstream ofAngle("outputs/variable_angle_thickness/angle.txt");

      for (double ll1 = var1min ; ll1 <= var1max+var1step*0.01 ; ll1+=var1step) {
        for (double ll2 = var2min ; ll2 <= var2max+var2step*0.01 ; ll2+=var2step) {
          if (var1 == 1) {
            angle = ll1*piDbl/180.;
            thicknesses.t[var1layerthick-1] = ll2;
            for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
            ml.setParameters(nbrMediums,Zi,mediums,lambda);
            ofAngle << ll1 << "\t";
            ofThickness << ll2*1e9 << "\t";
          }
          else {
            angle = ll2*piDbl/180.;
            thicknesses.t[var1layerthick-1] = ll1;
            for (i=1 ; i<nbrMediums-1 ; i++) Zi.t[i] = Zi.t[i-1]+thicknesses.t[i-1];
            ml.setParameters(nbrMediums,Zi,mediums,lambda);
            ofAngle << ll2 << "\t";
            ofThickness << ll1*1e9 << "\t";
          }
          if (angle<piDblby2) {
            ml.setKxy(ml.ki.t[0]*sin(angle));
            ampl.t[0] = 1.; ampl.t[1] = 0.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[0].re << "\t";
            ofSimr << aOut.t[0].im << "\t";
            ofSret << aOut.t[1].re << "\t";
            ofSimt << aOut.t[1].im << "\t";
            ofSR << norm(aOut.t[0]) << "\t";
            ofST << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[0].re << "\t";
            ofPimr << aOut.t[0].im << "\t";
            ofPret << aOut.t[1].re << "\t";
            ofPimt << aOut.t[1].im << "\t";
            ofPR << norm(aOut.t[0]) << "\t";
            ofPT << norm(aOut.t[1])*real(ml.calcKz(nbrMediums-1,true)/ml.calcKz(0,true)) << "\t";
          }
          else {
            ml.setKxy(ml.ki.t[nbrMediums-1]*sin(angle));
            ampl.t[0] = 0.; ampl.t[1] = 1.;
            aOut = ml.sourceExtOut(ampl,0);
            ofSrer << aOut.t[1].re << "\t";
            ofSimr << aOut.t[1].im << "\t";
            ofSret << aOut.t[0].re << "\t";
            ofSimt << aOut.t[0].im << "\t";
            ofSR << norm(aOut.t[1]) << "\t";
            ofST << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
            aOut = ml.sourceExtOut(ampl,1);
            ofPrer << aOut.t[1].re << "\t";
            ofPimr << aOut.t[1].im << "\t";
            ofPret << aOut.t[0].re << "\t";
            ofPimt << aOut.t[0].im << "\t";
            ofPR << norm(aOut.t[1]) << "\t";
            ofPT << norm(aOut.t[0])*real(ml.calcKz(0,false)/ml.calcKz(nbrMediums-1,false)) << "\t";
          }
          counting++;
        }
        cout << double(counting)/totCalc*100. << " %" << endl;
        ofAngle << endl;
        ofThickness << endl;
        ofSrer << endl;
        ofSimr << endl;
        ofSret << endl;
        ofSimt << endl;
        ofSR << endl;
        ofST << endl;
        ofPrer << endl;
        ofPimr << endl;
        ofPret << endl;
        ofPimt << endl;
        ofPR << endl;
        ofPT << endl;
      }
    }
  }
  cout << "\n Done!" << endl;
  return 0;
}
