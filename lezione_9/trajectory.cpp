#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

double Derivative_4Order(double (*func)(double), double x, double h);
double Derivative_FD(double (*func)(double), double x, double h);
double func(double x);
double Sec_Derivative(double (*func)(double), double x, double h);
double Sec_Derivative_fwd(double (*func)(double), double x, double h);

int main(){

    double alpha = 10.;
    int N = 100;
    double h = alpha/N;
    double x = 0.;
    

    ofstream fdata;         
    fdata.open("trajectory.dat");

    double traj = func(x);
    double vel = Derivative_FD(func, x, h);
    double acc = Sec_Derivative_fwd(func, x, h);

    while(x<alpha){
        fdata << x << " " << traj << " " << vel << " " << acc << endl;
        x += h;
        traj = func(x);
        vel = Derivative_FD(func, x, h);
        acc = Sec_Derivative_fwd(func, x, h);
    }

    fdata.close();

    return 0;
}

double func(double x){
    
    if(x == 0.) return 0.;
    else return 10.*x*x - x*x*x*(1 - exp(-(100.)/x));

}

double Derivative_4Order(double (*func)(double), double x, double h){

    double f_2prec, f_prec, f_2succ, f_succ;
    f_2prec = func(x - 2*h);
    f_prec = func(x - h);
    f_succ = func(x + h);
    f_2succ = func(x + 2*h);
    return (f_2prec - 8*f_prec + 8*f_succ - f_2succ)/(12*h);

}

double Derivative_FD(double (*func)(double), double x, double h){

    double f_succ, f_curr;
    f_succ = func(x + h);
    f_curr = func(x);
    return (f_succ - f_curr)/(2*h);

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
