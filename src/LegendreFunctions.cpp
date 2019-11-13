#include "LegendreFunctions.h"

LegendreFunctions::LegendreFunctions() : N(0) , pmn(NULL) , pimn(NULL) , taumn(NULL) {}

LegendreFunctions::LegendreFunctions(const int &NN,const double &xx) : N(NN) , x(xx) {
  if (NN >= 0) {
    pmn = new RVector[NN+1];
    pimn = new RVector[NN+1];
    taumn = new RVector[NN+1];
  }
  else return;
  if (x == 1.) calcElmts0();
  else if (x == -1.) calcElmtsPi();
  else calcElmts();
}

LegendreFunctions::~LegendreFunctions() {
  if (pmn != NULL) delete[] pmn;
  if (pimn != NULL) delete[] pimn;
  if (taumn != NULL) delete[] taumn;
}

void LegendreFunctions::assign(const int &NN, const double &xx) {
  N = NN;
  x = xx;
  if (pmn != NULL) delete [] pmn;
  if (pimn != NULL) delete [] pimn;
  if (taumn != NULL) delete [] taumn;
  if (NN >= 0) {
    pmn = new RVector[N+1];
    pimn = new RVector[N+1];
    taumn = new RVector[N+1];
  }
  else return;
  if (x == 1.) calcElmts0();
  else if (x == -1.) calcElmtsPi();
  else calcElmts();
}

void LegendreFunctions::calcElmts() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  double sint = sqrt(1.-x*x);
  if (N >= 0) {
    pmn[0].t[0] = sqrt(2.)*0.5;
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = sqrt(3.)*x*pmn[0].t[0];
    pimn[0].t[1] = pmn[0].t[1]/sint;
    taumn[0].t[1] = (x*pmn[0].t[1]-sqrt(3.)*pmn[0].t[0])/sint;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = sqrt((4.*double(i*i)-1.)/double(i*i))*(x*pmn[0].t[i-1]-sqrt((double(i-1)*double(i-1))/(4.*double(i-1)*double(i-1)-1.))*pmn[0].t[i-2]);
    pimn[0].t[i] = pmn[0].t[i]/sint;
    taumn[0].t[i] = (double(i)*x*pmn[0].t[i]-sqrt(double(2*i+1)*double(i*i)/double(2*i-1))*pmn[0].t[i-1])/sint;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = sqrt(double(2*i+1)/double(2*i))*sint*pmn[i-1].t[0];
    pimn[i].t[0] = pmn[i].t[0]/sint;
    taumn[i].t[0] = double(i)*x*pmn[i].t[0]/sint;
    if (N > i) {
      pmn[i].t[1] = sqrt(double(2*i+3))*x*pmn[i].t[0];
      pimn[i].t[1] = pmn[i].t[1]/sint;
      taumn[i].t[1] = (double(i+1)*x*pmn[i].t[1]-sqrt(double(2*(i+1)+1)*double((i+1)*(i+1)-i*i)/double(2*(i+1)-1))*pmn[i].t[0])/sint;
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = sqrt((4.*double((i+j)*(i+j))-1.)/double((i+j)*(i+j)-i*i))*(x*pmn[i].t[j-1]-sqrt((double(i+j-1)*double(i+j-1)-double(i*i))/(4.*double(i+j-1)*double(i+j-1)-1.))*pmn[i].t[j-2]);
      pimn[i].t[j] = pmn[i].t[j]/sint;
      taumn[i].t[j] = (double(i+j)*x*pmn[i].t[j]-sqrt(double(2*(i+j)+1)*double((i+j)*(i+j)-i*i)/double(2*(i+j)-1))*pmn[i].t[j-1])/sint;
    }
  }
}

void LegendreFunctions::calcElmts0() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  if (N >= 0) {
    pmn[0].t[0] = sqrt(0.5);
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = sqrt(1.5);
    pimn[0].t[1] = 0.;
    taumn[0].t[1] = 0.;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = sqrt(double(i+i+1)*0.5);
    pimn[0].t[i] = 0.;
    taumn[0].t[i] = 0.;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = 0.;
    if (i == 1) {
      pimn[i].t[0] = sqrt(3.)*0.5;
      taumn[i].t[0] = sqrt(3.)*0.5;
    }
    else {
      pimn[i].t[0] = 0.;
      taumn[i].t[0] = 0.;
    }
    if (N > i) {
      pmn[i].t[1] = 0.;
      if (i == 1) {
        pimn[i].t[1] = sqrt(15.)*0.5;
        taumn[i].t[1] = sqrt(15.)*0.5;
      }
      else {
        pimn[i].t[1] = 0.;
        taumn[i].t[1] = 0.;
      }
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = 0.;
      if (i == 1) {
        pimn[i].t[j] = sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
        taumn[i].t[j] = sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
      }
      else {
        pimn[i].t[j] = 0.;
        taumn[i].t[j] = 0.;
      }
    }
  }
}

void LegendreFunctions::calcElmtsPi() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  if (N >= 0) {
    pmn[0].t[0] = sqrt(0.5);
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = -sqrt(1.5);
    pimn[0].t[1] = 0.;
    taumn[0].t[1] = 0.;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = powm1(i)*sqrt(double(i+i+1)*0.5);
    pimn[0].t[i] = 0.;
    taumn[0].t[i] = 0.;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = 0.;
    if (i == 1) {
      pimn[i].t[0] = sqrt(3.)*0.5;
      taumn[i].t[0] = -sqrt(3.)*0.5;
    }
    else {
      pimn[i].t[0] = 0.;
      taumn[i].t[0] = 0.;
    }
    if (N > i) {
      pmn[i].t[1] = 0.;
      if (i == 1) {
        pimn[i].t[1] = -sqrt(15.)*0.5;
        taumn[i].t[1] = sqrt(15.)*0.5;
      }
      else {
        pimn[i].t[1] = 0.;
        taumn[i].t[1] = 0.;
      }
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = 0.;
      if (i == 1) {
        pimn[i].t[j] = powm1(j+2)*sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
        taumn[i].t[j] = -powm1(j+2)*sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
      }
      else {
        pimn[i].t[j] = 0.;
        taumn[i].t[j] = 0.;
      }
    }
  }
}

double LegendreFunctions::Pmn(const int &n, const int &m) {
  if (abs(m) > n) return 0.;
  if (m < 0) return powm1(-m)*pmn[-m].t[n+m];
  return pmn[m].t[n-m];
}
double LegendreFunctions::Pimn(const int &n, const int &m) {
  if (m < 0) return powm1(-m)*pimn[-m].t[n+m];
  return pimn[m].t[n-m];
}
double LegendreFunctions::Taumn(const int &n, const int &m) {
  if (m < 0) return powm1(-m)*taumn[-m].t[n+m];
  return taumn[m].t[n-m];
}




LegendreFunctionsC::LegendreFunctionsC() : N(0) , pmn(NULL) , pimn(NULL) , taumn(NULL) {}

LegendreFunctionsC::LegendreFunctionsC(const int &NN,const Complex &xx) : N(NN) , x(xx) {
  if (NN >= 0) {
    pmn = new Vector[NN+1];
    pimn = new Vector[NN+1];
    taumn = new Vector[NN+1];
  }
  else return;
  if (x == 1.) calcElmts0();
  else if (x == -1.) calcElmtsPi();
  else calcElmts();
}

LegendreFunctionsC::~LegendreFunctionsC() {
  if (pmn != NULL) delete[] pmn;
  if (pimn != NULL) delete[] pimn;
  if (taumn != NULL) delete[] taumn;
}

void LegendreFunctionsC::assign(const int &NN, const Complex &xx) {
  N = NN;
  x = xx;
  if (pmn != NULL) delete [] pmn;
  if (pimn != NULL) delete [] pimn;
  if (taumn != NULL) delete [] taumn;
  if (NN >= 0) {
    pmn = new Vector[N+1];
    pimn = new Vector[N+1];
    taumn = new Vector[N+1];
  }
  else return;
  if (x == 1.) calcElmts0();
  else if (x == -1.) calcElmtsPi();
  else calcElmts();
}

void LegendreFunctionsC::calcElmts() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  Complex sint = sqrt(1.-x*x);
  if (N >= 0) {
    pmn[0].t[0] = sqrt(2.)*0.5;
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = sqrt(3.)*x*pmn[0].t[0];
    pimn[0].t[1] = pmn[0].t[1]/sint;
    taumn[0].t[1] = (x*pmn[0].t[1]-sqrt(3.)*pmn[0].t[0])/sint;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = sqrt((4.*double(i*i)-1.)/double(i*i))*(x*pmn[0].t[i-1]-sqrt((double(i-1)*double(i-1))/(4.*double(i-1)*double(i-1)-1.))*pmn[0].t[i-2]);
    pimn[0].t[i] = pmn[0].t[i]/sint;
    taumn[0].t[i] = (double(i)*x*pmn[0].t[i]-sqrt(double(2*i+1)*double(i*i)/double(2*i-1))*pmn[0].t[i-1])/sint;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = sqrt(double(2*i+1)/double(2*i))*sint*pmn[i-1].t[0];
    pimn[i].t[0] = pmn[i].t[0]/sint;
    taumn[i].t[0] = double(i)*x*pmn[i].t[0]/sint;
    if (N > i) {
      pmn[i].t[1] = sqrt(double(2*i+3))*x*pmn[i].t[0];
      pimn[i].t[1] = pmn[i].t[1]/sint;
      taumn[i].t[1] = (double(i+1)*x*pmn[i].t[1]-sqrt(double(2*(i+1)+1)*double((i+1)*(i+1)-i*i)/double(2*(i+1)-1))*pmn[i].t[0])/sint;
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = sqrt((4.*double((i+j)*(i+j))-1.)/double((i+j)*(i+j)-i*i))*(x*pmn[i].t[j-1]-sqrt((double(i+j-1)*double(i+j-1)-double(i*i))/(4.*double(i+j-1)*double(i+j-1)-1.))*pmn[i].t[j-2]);
      pimn[i].t[j] = pmn[i].t[j]/sint;
      taumn[i].t[j] = (double(i+j)*x*pmn[i].t[j]-sqrt(double(2*(i+j)+1)*double((i+j)*(i+j)-i*i)/double(2*(i+j)-1))*pmn[i].t[j-1])/sint;
    }
  }
}

void LegendreFunctionsC::calcElmts0() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  if (N >= 0) {
    pmn[0].t[0] = sqrt(0.5);
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = sqrt(1.5);
    pimn[0].t[1] = 0.;
    taumn[0].t[1] = 0.;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = sqrt(double(i+i+1)*0.5);
    pimn[0].t[i] = 0.;
    taumn[0].t[i] = 0.;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = 0.;
    if (i == 1) {
      pimn[i].t[0] = sqrt(3.)*0.5;
      taumn[i].t[0] = sqrt(3.)*0.5;
    }
    else {
      pimn[i].t[0] = 0.;
      taumn[i].t[0] = 0.;
    }
    if (N > i) {
      pmn[i].t[1] = 0.;
      if (i == 1) {
        pimn[i].t[1] = sqrt(15.)*0.5;
        taumn[i].t[1] = sqrt(15.)*0.5;
      }
      else {
        pimn[i].t[1] = 0.;
        taumn[i].t[1] = 0.;
      }
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = 0.;
      if (i == 1) {
        pimn[i].t[j] = sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
        taumn[i].t[j] = sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
      }
      else {
        pimn[i].t[j] = 0.;
        taumn[i].t[j] = 0.;
      }
    }
  }
}

void LegendreFunctionsC::calcElmtsPi() {
  int i,j;
  for (int i=0 ; i<=N ; i++) {
    pmn[i].resize(N+1-i);
    pimn[i].resize(N+1-i);
    taumn[i].resize(N+1-i);
  }
  if (N >= 0) {
    pmn[0].t[0] = sqrt(0.5);
    pimn[0].t[0] = 0.;
    taumn[0].t[0] = 0.;
  }
  if (N > 0) {
    pmn[0].t[1] = -sqrt(1.5);
    pimn[0].t[1] = 0.;
    taumn[0].t[1] = 0.;
  }
  for (i=2 ; i<=N ; i++) {
    pmn[0].t[i] = powm1(i)*sqrt(double(i+i+1)*0.5);
    pimn[0].t[i] = 0.;
    taumn[0].t[i] = 0.;
  }
  for (i=1 ; i<=N ; i++) {
    pmn[i].t[0] = 0.;
    if (i == 1) {
      pimn[i].t[0] = sqrt(3.)*0.5;
      taumn[i].t[0] = -sqrt(3.)*0.5;
    }
    else {
      pimn[i].t[0] = 0.;
      taumn[i].t[0] = 0.;
    }
    if (N > i) {
      pmn[i].t[1] = 0.;
      if (i == 1) {
        pimn[i].t[1] = -sqrt(15.)*0.5;
        taumn[i].t[1] = sqrt(15.)*0.5;
      }
      else {
        pimn[i].t[1] = 0.;
        taumn[i].t[1] = 0.;
      }
    }
    for (j=2 ; j<=N-i ; j++) {
      pmn[i].t[j] = 0.;
      if (i == 1) {
        pimn[i].t[j] = powm1(j+2)*sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
        taumn[i].t[j] = -powm1(j+2)*sqrt(double((j+1)*(j+1)+(j+1))*double((j+1)+(j+1)+1)*0.5)*0.5;
      }
      else {
        pimn[i].t[j] = 0.;
        taumn[i].t[j] = 0.;
      }
    }
  }
}

Complex LegendreFunctionsC::Pmn(const int &n, const int &m) {
  if (abs(m) > n) return 0.;
  if (m < 0) return powm1(-m)*pmn[-m].t[n+m];
  return pmn[m].t[n-m];
}
Complex LegendreFunctionsC::Pimn(const int &n, const int &m) {
  if (m < 0) return powm1(-m)*pimn[-m].t[n+m];

  return pimn[m].t[n-m];
}
Complex LegendreFunctionsC::Taumn(const int &n, const int &m) {
  if (m < 0) return powm1(-m)*taumn[-m].t[n+m];
  return taumn[m].t[n-m];
}
