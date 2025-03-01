#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
double Secant(double (*func)(double), double a, double b, double tol, int& ntry);
double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry);
double FalsePos(double (*func)(double), double a, double b, double tol, int& ntry);
double Horner_pol(const double, const double[], const int, double&);
void Bracketing(double (*func)(double), double a, double b, int N, double XL[], double XR[], int& k);