#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

double func(double x);
double dfunc(double x);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
double Secant(double (*func)(double), double a, double b, double tol, int& ntry);
double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry);
double FalsePos(double (*func)(double), double a, double b, double tol, int& ntry);

int main(){

    int ntry = 0;
    cout << setiosflags(ios::scientific) << setprecision(8);
    cout << "Bisection: \t" <<Bisection(func, -1., 1., 1.e-7, ntry) << "  ntry: \t" << ntry << endl;
    ntry = 0;
    cout << "FalsePos: \t" << FalsePos(func, -1., 1., 1.e-7, ntry) << "  ntry: \t" << ntry <<endl;
    cout << "===================================" << endl;
    cout << "Secant: \t" << Secant(func, -1., 1., 1.e-7, ntry) << "  ntry: \t" << ntry << endl;
    ntry = 0;
    cout << "Newton: \t" << Newton(func, dfunc, -1., 1., 1.e-7, ntry) << "  ntry: \t" << ntry << endl;


    return 0;
}

double func(double x){
    return exp(-x) - x;
}

double dfunc(double x){
    return -exp(-x) - 1;
}

double Bisection(double (*func)(double), double a, double b, double tol, int& ntry){

    double xm, f_xm;
    double sup_lim, inf_lim;
    double f_inf, f_sup;

    sup_lim = b;
    inf_lim = a;
    xm = (a + b)/2.;
    f_xm = func(xm);
    f_inf = func(a);
    f_sup = func(b);
    ntry = 1;

    if(f_sup == 0.) return sup_lim;
    else if(f_inf == 0.) return inf_lim;

    while(fabs(sup_lim - inf_lim) > tol){

        if(f_xm == 0.) break;
        else if(f_inf * f_xm < 0){
            sup_lim = xm;
            xm = (sup_lim + inf_lim) / 2.;
            f_xm = func(xm);
        }
        else if(f_inf * f_xm > 0){
            inf_lim = xm;
            f_inf = f_xm;
            xm = (sup_lim + inf_lim) / 2.;
            f_xm = func(xm); 
        }
        ntry ++;
    }

    return xm;
}


double FalsePos(double (*func)(double), double a, double b, double tol, int& ntry){

    double x1, x2;
    double f_inf, f_sup;
    double inf_lim, sup_lim;
    double m, c;

    inf_lim = a;
    sup_lim = b;
    f_inf = func(a);
    f_sup = func(b);
    m = (f_sup - f_inf) / (sup_lim - inf_lim);
    c = f_inf - m*inf_lim;
    x1 = a;
    x2 = - c/m;
    ntry = 1;

    if(f_sup == 0.) return sup_lim;
    else if(f_inf == 0.) return inf_lim;

    while(fabs(x2 - x1) > tol){
        x1 = x2;
        inf_lim = x1;
        f_inf = func(x1);
        m = (f_sup - f_inf) / (sup_lim - inf_lim);
        c = f_inf - m*inf_lim;
        x2 = - c/m;
        ntry ++;
    }

    return x2;

}

double Secant(double (*func)(double), double a, double b, double tol, int& ntry){

    double sup_lim, inf_lim;
    double f_inf, f_sup;
    double dx;

    inf_lim = a;
    sup_lim = b;
    dx = sup_lim - inf_lim;
    f_inf = func(inf_lim);
    f_sup = func(sup_lim);
    ntry = 1;

    while(fabs(dx) > tol){
        dx = f_sup * (sup_lim - inf_lim)/(f_sup - f_inf);    
        inf_lim = sup_lim;
        f_inf = f_sup;
        sup_lim = sup_lim - dx;
        f_sup = func(sup_lim);
        
        ntry ++;
        }

    return sup_lim;
}


double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry){

    double sup_lim, inf_lim;
    double f, f_der;
    double xc, x_guess, dx;

    inf_lim = a;
    sup_lim = b;
    xc = (sup_lim + inf_lim)/2.;
    f = func(xc);
    f_der = dfunc(xc);
    x_guess = xc - f/f_der;
    dx = x_guess - xc;
    ntry = 1;

    while(fabs(dx) > tol){ 
        xc = x_guess;
        f = func(x_guess);
        f_der = dfunc(x_guess);
        dx = - f/f_der;
        x_guess = xc + dx;
        
        ntry ++;
        }

    return x_guess;
}

