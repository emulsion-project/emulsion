#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <cmath>
using namespace std;

static const double piDbl = 3.14159265358979323846; // pi
static const double piDblby2 = piDbl / 2.; // pi/2
static const double piDbl2 = 2. * piDbl; //2*pi
static const double me = 9.1093826e-31; // electron mass (kg)
static const double e = 1.602176487e-19; // elementary charge (C)
static const double eps0 = 8.854187e-12; // free space dielectric permittivity (As/Vm)
static const double mu0 = 4. * piDbl * 1.0e-7; // free space magnetic permeability
static const double c = 299792458;  // free space light speed (m/s)
static const double h = 6.62606896e-34;  // Planck constant (J.s)
static const double hbar = 1.054571628e-34; // reduced Plank constant h/2pi (J.s)
static const double euler = 0.5772156649015328606;

int ii(const int &n, const int &m);
double powm1(const int &n);
void swap(int &a, int &b);
void swap(double &a, double &b);
double max(const double &a, const double &b);
int max(const int &a, const int &b);
double fact(const int &n); // compute n!
double dblfact(const int &n); // compute n!!
double binomialCoef(const int &n, const int &l);

#endif // UTILS_H_INCLUDED
