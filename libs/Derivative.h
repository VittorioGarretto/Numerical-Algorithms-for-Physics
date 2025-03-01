#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double Derivative_FD(double (*func)(double), double x, double h);
double Derivative_BD(double (*func)(double), double x, double h);
double Derivative_CD(double (*func)(double), double x, double h); 
double Derivative_4Order(double (*func)(double), double x, double h); 
double Sec_Derivative(double (*func)(double), double x, double h);
double Sec_Derivative_fwd(double (*func)(double), double x, double h);