#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void pVerlet(double *Y, double *v, void (*Acceleration)(double *, double *), double dt, int neq);
