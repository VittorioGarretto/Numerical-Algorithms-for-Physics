#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

double Derivative_FD(double (*func)(double), double x, double h);
double Derivative_BD(double (*func)(double), double x, double h);
double Derivative_CD(double (*func)(double), double x, double h); 
double Derivative_4Order(double (*func)(double), double x, double h); 
double func(double x);
double dfunc(double x);

int main(){

    double x = 1.;
    double h = 0.5;
    double tol = 1.e-12;
    double f_der_ex = dfunc(x); 
    double f_der;

    ofstream fdata;         
    fdata.open("FD.dat");
    f_der = Derivative_FD(func, x, h);
    while(1./h < 1000){
        fdata << 1/h << " " << fabs(f_der - f_der_ex) << endl;
        h /= 2.;
        f_der = Derivative_FD(func, x, h);
    }
    fdata.close();

    h = 0.5;
    fdata.open("BD.dat");
    f_der = Derivative_BD(func, x, h);
    while(1./h < 1000){
        fdata << 1/h << " " << fabs(f_der - f_der_ex) << endl;
        h /= 2.;
        f_der = Derivative_BD(func, x, h);
    }
    fdata.close();

    h = 0.5;
    fdata.open("CD.dat");
    f_der = Derivative_CD(func, x, h);
    while(1./h < 1000){
        fdata << 1/h << " " << fabs(f_der - f_der_ex) << endl;
        h /= 2.;
        f_der = Derivative_CD(func, x, h);
    }
    fdata.close();

    h = 0.5;
    fdata.open("4Order.dat");
    f_der = Derivative_4Order(func, x, h);
    while(1./h < 1000){
        fdata << 1/h << " " << fabs(f_der - f_der_ex) << endl;
        h /= 2.;
        f_der = Derivative_4Order(func, x, h);
    }
    fdata.close();

    return 0;
}

double func(double x){
    return sin(x);
}

double dfunc(double x){
    return cos(x);
}

double Derivative_FD(double (*func)(double), double x, double h){

    double f_succ, f_curr;
    f_succ = func(x + h);
    f_curr = func(x);
    return (f_succ - f_curr)/(h);

}

double Derivative_BD(double (*func)(double), double x, double h){

    double f_prec, f_curr;
    f_prec = func(x - h);
    f_curr = func(x);
    return (f_curr - f_prec)/(h);

}

double Derivative_CD(double (*func)(double), double x, double h){

    double f_prec, f_succ;
    f_prec = func(x - h);
    f_succ = func(x + h);
    return (f_succ - f_prec)/(2*h);

}

double Derivative_4Order(double (*func)(double), double x, double h){

    double f_2prec, f_prec, f_2succ, f_succ;
    f_2prec = func(x - 2*h);
    f_prec = func(x - h);
    f_succ = func(x + h);
    f_2succ = func(x + 2*h);
    return (f_2prec - 8*f_prec + 8*f_succ - f_2succ)/(12*h);

}
