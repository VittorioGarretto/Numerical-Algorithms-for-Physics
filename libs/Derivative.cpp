#include "Derivative.h"

double Derivative_FD(double (*func)(double), double x, double h){

    double f_succ, f_curr;
    f_succ = func(x + h);
    f_curr = func(x);
    return (f_succ - f_curr)/h;

}

double Derivative_BD(double (*func)(double), double x, double h){

    double f_prec, f_curr;
    f_prec = func(x - h);
    f_curr = func(x);
    return (f_curr - f_prec)/h;

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

double Sec_Derivative_fwd(double (*func)(double), double x, double h) {

    double f_curr = func(x);
    double f_succ = func(x + h);
    double f_2succ = func(x + 2 * h);
    return (f_2succ - 2*f_succ + f_curr) / (h*h);

}

double Sec_Derivative(double (*func)(double), double x, double h){

    double f_curr, f_succ, f_prec;
    f_curr = func(x);
    f_succ = func(x + h);
    f_prec = func(x - h);
    return (f_succ - 2*f_curr + f_prec)/(h*h);

}