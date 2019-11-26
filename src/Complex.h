#ifndef COMPLEX_H_INCLUDED
#define COMPLEX_H_INCLUDED

#include <ostream>
#include "Utils.h"

class Complex {
public:
	Complex(void);
	Complex(const double &x);
	Complex(const Complex &z);
	Complex(const double &a, const double &b);
	double real() const;
	double imag() const;
	Complex& operator=(const Complex &z);
	Complex& operator=(const double &x);
	Complex& operator+=(const Complex &z);
	Complex& operator+=(const double &x);
	Complex& operator-=(const Complex &z);
	Complex& operator-=(const double &x);
	Complex& operator*=(const Complex &z);
	Complex& operator*=(const double &x);
	Complex& operator/=(const Complex &z);
	Complex& operator/=(const double &x);
	double re, im;
};

/*********************************
 *        Related operators      *
 *********************************/

Complex operator+(const Complex &z1, const Complex &z2);
Complex operator+(const Complex &z, const double &x);
Complex operator+(const double &x, const Complex &z);
Complex operator+(const Complex &z);

/************************/

Complex operator-(const Complex &z1, const Complex &z2);
Complex operator-(const Complex &z, const double &x);
Complex operator-(const double &x, const Complex &z);
Complex operator-(const Complex &z);

/************************/

Complex operator*(const Complex &z1, const Complex &z2);
Complex operator*(const Complex &z, const double &x);
Complex operator*(const double &x, const Complex &z);

/************************/

Complex operator/(const Complex &z1, const Complex &z2);
Complex operator/(const Complex &z, const double &x);
Complex operator/(const double &x, const Complex &z);

/************************/

bool operator==(const Complex &z1, const Complex &z2);
bool operator==(const Complex &z, const double &x);
bool operator==(const double &x, const Complex &z);

/************************/

bool operator!=(const Complex &z1, const Complex &z2);
bool operator!=(const Complex &z, const double &x);
bool operator!=(const double &x, const Complex &z);

/************************/

bool operator<(const Complex &z1, const Complex &z2);
bool operator<=(const Complex &z1, const Complex &z2);

/************************/

bool operator>(const Complex &z1, const Complex &z2);
bool operator>=(const Complex &z1, const Complex &z2);

/*********************************
 *        Definition of i,j      *
 *********************************/

static const Complex i_(0., 1.);
static const Complex j_(0., 1.);

/*********************************
 *      Related functions        *
 *********************************/

Complex cart(const double &a, const double &b);
double real(const Complex &z);
double imag(const Complex &z);
double abs(const Complex &z);
double norm(const Complex &z);
double phase(const Complex &z);
Complex conj(const Complex &z);
Complex polar(const double &r, const double &phi);

/************************/

Complex sqrt(const Complex &z);
Complex exp(const Complex &z);
Complex log(const Complex &z);
Complex log10(const Complex &z);
Complex pow(const Complex &z, const double &x);
Complex pow(const double &x, const Complex &z);
Complex pow(const Complex &z1, const Complex &z2);
Complex powi(const int &n);

/************************/

Complex sin(const Complex &z);
Complex cos(const Complex &z);
Complex tan(const Complex &z);
Complex cotan(const Complex &z);
Complex asin(const Complex &z);
Complex acos(const Complex &z);
Complex atan(const Complex &z);
Complex cosh(const Complex &z);
Complex sinh(const Complex &z);
Complex tanh(const Complex &z);
Complex acosh(const Complex &z);
Complex asinh(const Complex &z);
Complex atanh(const Complex &z);

/************************/

void swap(Complex &a, Complex &b);

/*********************************
 *      Stream operators         *
 *********************************/

std::ostream& operator<<(std::ostream &out, const Complex &z);

#endif // COMPLEX_H_INCLUDED
