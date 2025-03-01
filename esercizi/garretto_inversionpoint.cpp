// Name: Vittorio, Garretto
// Date: 14 nov 2024
// 
// Code output:
// **********************************************************
// tlo = 6.000000e+00; thi = 7.000000e+00; nint= 2.000000e+00; Func eval = 10
// tlo = 6.000000e+00; thi = 7.765933e+00; nint= 4.000000e+00; Func eval = 30
// tlo = 6.000000e+00; thi = 7.765933e+00; nint= 4.000000e+00; Func eval = 50
// tlo = 6.000000e+00; thi = 7.858601e+00; nint= 4.000000e+00; Func eval = 70
// tlo = 6.000000e+00; thi = 7.858601e+00; nint= 4.000000e+00; Func eval = 90
// tlo = 6.000000e+00; thi = 7.859718e+00; nint= 4.000000e+00; Func eval = 110
// 
// inversion time (zero) = 7.859718e+00
// **********************************************************


#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

double Gauss (double (*F)(double), double a, double b, int N, int Ng);
double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry);
double Velocity(double t);
double Acceleration(double t);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(6);
    double tol = 1.e-6;
    int ntry = 0;

    double v0 = -3.937816e-02; // ottimizzazione trovata calcolando Velocity(6.) e assumendo tlo = 6.

    cout << Newton(Velocity, Acceleration, 6., 8., tol, ntry) << endl;

    return 0;
}


double Velocity(double t){

    static int nfv = 0;
    double tlo = 6.; // ottimizzazione
    double thi = t;
    double v0 = -3.937816e-02; // ottimizzazione (nuova v(0))
    int nint = ceil(fabs(tlo - thi)*2);
    double v_t;
    int Ng = 5;

    v_t = Gauss(Acceleration, tlo, thi, nint, Ng);
    nfv += Ng*nint;

    cout << "tlo = " << tlo << "; thi = "<< t << "; nint = "<< 
    nint << "; Func eval = " << nfv << endl;
    
    return v_t + v0;

}

double Acceleration(double t){

    return (1. - exp(-t))/(sqrt(1. + t*t*t*t));
    
}

double Gauss(double (*F)(double), double a, double b, int N, int Ng){

    double w[8];
    double x[8];

    if (Ng == 2) {
        x[0] =  sqrt(1.0 / 3.0);     w[0] = 1;
        x[1] = -sqrt(1.0 / 3.0);     w[1] = 1;
    } 
    else if (Ng == 3) {
        x[0] =  0;						 w[0] = 8.0 / 9.0;
        x[1] =  sqrt(3.0 / 5.0);		 w[1] = 5.0 / 9.0;
        x[2] = -sqrt(3.0 / 5.0);		 w[2] = 5.0 / 9.0;
    }
    else if (Ng == 4){
        x[0] =  sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[0] = (18 + sqrt(30)) / 36;
        x[1] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[1] = (18 + sqrt(30)) / 36;
        x[2] =  sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[2] = (18 - sqrt(30)) / 36;
        x[3] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[3] = (18 - sqrt(30)) / 36;
    }
    else if (Ng == 5){
        x[0] =  0;												w[0] = 128.0 / 225.0;
		x[1] =  1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		w[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0; 
		x[2] = -1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		w[2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		x[3] =  1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		w[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		x[4] = -1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		w[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
    }

    double sum = 0.;
    double xc;  // centro di ciascun sottointervallo
    double h = (b-a)/N;
    double sumk;

    for(int i=0;i<N;i++){  // Ciclo su ciascuno dei N sottointervalli
        xc = ((a + i*h) + (a + (i+1)*h))/2.;
        sumk = 0.;
        for(int k=0;k<Ng;k++) {
            sumk += (w[k] * F(h*0.5*x[k] + xc)); // Regola di gauss sui sottointervalli
        } 
        sum += sumk;
    }

    return sum*h*0.5;
    
}

double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry) {

    double sup_lim = b, inf_lim = a;   
    double xc = (sup_lim + inf_lim) / 2.0;  
    double f = func(xc);
    double f_der = dfunc(xc);

    double x_guess = xc - f/f_der;   // Prima stima della radice
    double dx = x_guess - xc;
    ntry = 1;  

    while (fabs(func(x_guess)) > tol) { 
        xc = x_guess;                  
        f = func(xc);         
        f_der = dfunc(xc); 
        dx = -f / f_der;               
        x_guess = xc + dx;
        ntry ++;
    }

    return x_guess;  
}