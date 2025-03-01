#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double Rectangular(double (*F)(double), double a, double b, int N);
double Trapezoidal(double (*F)(double), double a, double b, int N);
double Simpson(double (*F)(double), double a, double b, int N);
double Gauss (double (*F)(double), double a, double b, int N, int Ng);
double Gaussian2DQuad (double (*Func)(double, double), double xa, double xb, double ya, double yb, int N, int Ng);