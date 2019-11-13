#include "Complex.h"
#include <iostream>
//using namespace std;

/*********************************
 *     Class constructors        *
 *********************************/

Complex::Complex() : re(0.) , im(0.) {}
Complex::Complex(const double &x) : re(x) , im(0.) {}
Complex::Complex(const double &x, const double &y) : re(x) , im(y) {}
Complex::Complex(const Complex &z) { re = z.re; im = z.im; }

/*********************************
 *     Class member functions    *
 *********************************/

double Complex::real() const { return re; }
double Complex::imag() const { return im; }

/*********************************
 *     Class member operators    *
 *********************************/

Complex& Complex::operator=(const Complex &z) {
  if(&z!=this) { re = z.re; im = z.im; }
  return *this;
}

Complex& Complex::operator=(const double &x) {
  this->re = x;
  this->im = 0.;
  return *this;
}

Complex& Complex::operator+=(const Complex &z) {
  re += z.re;
  im += z.im;
  return *this;
}

Complex& Complex::operator+=(const double &x) {
  re += x;
  return *this;
}

Complex& Complex::operator-=(const Complex &z) {
  re -= z.re;
  im -= z.im;
  return *this;
}

Complex& Complex::operator-=(const double &x) {
  re -= x;
  return *this;
}

Complex& Complex::operator*=(const Complex &z) {
  double r = re;
  re = re*z.re - im*z.im;
  im = r*z.im + im*z.re;
  return *this;
}

Complex& Complex::operator*=(const double &x) {
  re *= x;
  im *= x;
  return *this;
}

Complex& Complex::operator/=(const Complex &z) {
  double r = re , d = z.re*z.re + z.im*z.im;
  re = (re*z.re + im*z.im) / d;
  im = (im*z.re - r*z.im) / d;
  return *this;
}

Complex& Complex::operator/=(const double &x) {
  re /= x;
  im /= x;
  return *this;
}

/*********************************
 *        Related operators      *
 *********************************/

Complex operator+(const Complex &z1, const Complex &z2) { return Complex(z1.re+z2.re, z1.im+z2.im); }
Complex operator+(const Complex &z, const double &x) { return Complex(z.re+x,z.im); }
Complex operator+(const double &x, const Complex &z) { return Complex(z.re+x,z.im); }
Complex operator+(const Complex &z) { return z; }

/************************/

Complex operator-(const Complex &z1, const Complex &z2) { return Complex(z1.re-z2.re,z1.im-z2.im); }
Complex operator-(const Complex &z, const double &x) { return Complex(z.re-x,z.im); }
Complex operator-(const double &x, const Complex &z) { return Complex(x-z.re,-z.im); }
Complex operator-(const Complex &z) { return Complex(-z.re,-z.im); }

/************************/

Complex operator*(const Complex &z1, const Complex &z2) {
  return Complex(z1.re*z2.re-z1.im*z2.im,z1.re*z2.im+z1.im*z2.re);
}
Complex operator*(const Complex &z, const double &x) { return Complex(z.re*x,z.im*x); }
Complex operator*(const double &x, const Complex &z) { return Complex(z.re*x,z.im*x); }

/************************/

Complex operator/(const Complex &z1, const Complex &z2) {
  double m = z2.re*z2.re+z2.im*z2.im;
  return Complex((z1.re*z2.re+z1.im*z2.im)/(m),(z1.im*z2.re-z1.re*z2.im)/(m));
}

Complex operator/(const Complex &z, const double &x) { return Complex(z.re/x,z.im/x); }

Complex operator/(const double &x, const Complex &z) {
  double norm = z.re*z.re + z.im*z.im;
  return Complex(z.re*x/norm,-z.im*x/norm);
}

/************************/

bool operator==(const Complex &z1, const Complex &z2) {
  if (z1.re==z2.re && z1.im==z2.im) return true;
  else return false;
}

bool operator==(const Complex &z, const double &x) {
  if ((z.re == x) && (z.im == 0.)) return true;
  else return false;
}

bool operator==(const double &x, const Complex &z) {
  if ((z.re == x) && (z.im == 0.)) return true;
  else return false;
}

/************************/

bool operator!=(const Complex &z1, const Complex &z2) {
  if ((z1.re == z2.re) && (z1.im == z2.im)) return false;
  else return true;
}

bool operator!=(const Complex &z, const double &x) {
  if ((z.re == x) && (z.im == 0)) return false;
  else return true;
}

bool operator!=(const double &x, const Complex &z) {
  if ((z.re==x) && (z.im==0.0)) return false;
  else return true;
}

/************************/

bool operator<(const Complex &z1, const Complex &z2) {
  double abs1 = sqrt(z1.re*z1.re + z1.im*z1.im),
    abs2 = sqrt(z2.re*z2.re + z2.im*z2.im);
  if (abs1 < abs2) return true;
  else return false;
}

bool operator<=(const Complex &z1, const Complex &z2) {
  double abs1 = sqrt(z1.re*z1.re + z1.im*z1.im),
    abs2 = sqrt(z2.re*z2.re + z2.im*z2.im);
  if (abs1 <= abs2) return true;
  else return false;
}

/************************/

bool operator>(const Complex &z1, const Complex &z2) {
  double abs1 = sqrt(z1.re*z1.re + z1.im*z1.im),
    abs2 = sqrt(z2.re*z2.re + z2.im*z2.im);
  if (abs1 > abs2) return true;
  else return false;
}

bool operator>=(const Complex &z1, const Complex &z2) {
  double abs1 = sqrt(z1.re*z1.re + z1.im*z1.im),
    abs2 = sqrt(z2.re*z2.re + z2.im*z2.im);
  if (abs1 > abs2) return true;
  else return false;
}

/*********************************
 *      Related functions        *
 *********************************/

Complex cart(const double &a, const double &b) { return Complex(a,b); }
double real(const Complex &z) { return z.re; }
double imag(const Complex &z) { return z.im; }
double abs(const Complex &z) { return sqrt(z.re*z.re + z.im*z.im); }
double norm(const Complex &z) { return (z.re*z.re + z.im*z.im); }
double phase(const Complex &z) { return atan2(z.im,z.re); }
Complex conj(const Complex &z) { return Complex(z.re,-z.im); }
Complex polar(const double &r, const double &phi) { return Complex(r*cos(phi),r*sin(phi)); }

/************************/

Complex sqrt(const Complex &z) { return polar(sqrt(abs(z)),phase(z)/2.); }
Complex exp(const Complex &z) { return cart(exp(z.re)*std::cos(z.im),exp(z.re)*std::sin(z.im)); }
Complex log(const Complex &z) { return cart(log(abs(z)),phase(z)); }
Complex log10(const Complex &z) { return log(z)/log(10.); }
Complex pow(const Complex &z, const double &x) { return polar(pow(abs(z),x),phase(z)*x); }
Complex pow(const double &x, const Complex &z) { return polar(pow(abs(z),x),phase(z)*x); }
Complex pow(const Complex &z1, const Complex &z2) { return exp(z2*log(z1)); }
Complex powi(const int &n) {
    if ((n%2)==0) return Complex(real(pow(i_,n)),0.);
    return Complex(0.,imag(pow(i_,n)));
}

/************************/

Complex sin(const Complex &z) { return cart(sin(z.re)*cosh(z.im),cos(z.re)*sinh(z.im)); }
Complex cos(const Complex &z) { return cart(cos(z.re)*cosh(z.im),-sin(z.re)*sinh(z.im)); }
Complex tan(const Complex &z) { return sin(z)/cos(z); }
Complex cotan(const Complex &z) { return cos(z)/sin(z); }
Complex asin(const Complex &z) { return -i_*log(i_*z+sqrt(1.-z*z)); }
Complex acos(const Complex &z) { return -i_*log(z+sqrt(z*z-1.)); }
Complex atan(const Complex &z) { return i_*log((1.-i_*z)/(1.+i_*z))/2.; }
Complex cosh(const Complex &z) { return cart(cos(z.im)*cosh(z.re),sin(z.im)*sinh(z.re)); }
Complex sinh(const Complex &z) { return cart(sinh(z.re)*cos(z.im),cosh(z.re)*sin(z.im)); }
Complex tanh(const Complex &z) { return sinh(z)/cosh(z); }
Complex acosh(const Complex &z) { return log(z+sqrt(z*z-1.)); }
Complex asinh(const Complex &z) { return log(z+sqrt(z*z+1.)); }
Complex atanh(const Complex &z) { return log((1.+z)/(1.-z))*0.5; }

/************************/

void swap(Complex &a,Complex &b) {
  Complex c = a;
  a=b;
  b=c;
}

/*********************************
 *      Stream operators         *
 *********************************/

std::ostream &operator<<(std::ostream &out, const Complex &z) {
  if (z.im == 0.) out << z.re;
  else out << '(' << z.re << ',' << z.im << ')';
  return out;
}
