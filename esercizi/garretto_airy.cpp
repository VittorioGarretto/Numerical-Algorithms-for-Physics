// Name: Vittorio, Garretto
// Date: 12 dec 2024
// 
// Code output:
// ********************************************************** 
// y' (xl) = 2.475674e+01
// **********************************************************

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void TridiagonalSolver(int n, double a[], double b[], double c[], double r[], double x[]);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
void dYdt(double x, double *Y, double *R);
double Residual (double s);


int main(){

    /////////////// Runge-Kutta /////////////////////

    double Y[2];
    int neq = 2;
    int nsteps = 800;
    double xl = -10.;
    double xr = 10.;
    double dx = fabs(xr - xl)/double(nsteps);

    double s1 = 10.;
    double s2 = 40.;
    int ntry = 0;
    double tol = 1.e-8;

    cout << setiosflags(ios::scientific) << setprecision(6);

    double s_true = Bisection(Residual, s1, s2, tol, ntry);
    cout << "Bisection -> s = " << s_true << endl;

    ofstream fdata;         
    fdata.open("airy.dat");
    Y[0] = 1.;
    Y[1] = s_true;
    fdata << xl << " " << Y[0] << endl;
    for(int i=0;i<nsteps;i++){
        RK4Step(xl, Y, dYdt, dx, neq);
        xl += dx;
        fdata << xl << " " << Y[0] << endl;
    }
    fdata << endl << endl;

    /////////////// Tridiagonal-Solver /////////////////////

    int n = nsteps + 1; 
    xl = -10.;  
    dx = fabs(xr - xl)/double(n-1);
    double x0 = 1.;
    double x1 = 0.;
    double a[n], b[n], c[n], x[n], r[n];
    x[0] = x0;
    x[n-1] = x1;

    for(int i=0;i<n;i++) a[i] = 1.;
    a[0] = a[1] = NAN;

    for(int i=0;i<n;i++) b[i] = -2. - (xl + i*dx)*dx*dx;
    b[0] = b[n-1] = NAN;

    for(int i=0;i<n;i++) c[i] = 1.;
    c[n-1] = c[n-2] = NAN;

    for(int i=0;i<n;i++) r[i] = 0.;
    r[1] = -x0;
    r[n-2] = -x1;
    r[0] = r[n-1] = NAN;

    TridiagonalSolver(n-2, a+1, b+1, c+1, r+1, x+1);

    for (int i=0;i<n;i++) {
        fdata << xl + i*dx << " " << x[i] << endl;
    }
    fdata.close();

    return 0;
}


void dYdt(double x, double *Y, double *R){

    R[0] = Y[1];
    R[1] = x*Y[0];

}

double Residual(double s){

    double Y[2];
    int neq = 2;
    int nsteps = 800;
    double xl = -10.;
    double xr = 10.;
    double dx = fabs(xr - xl)/double(nsteps);

    Y[0] = 1.;
    Y[1] = s;

    for(int i=0;i<nsteps;i++){
        RK4Step(xl, Y, dYdt, dx, neq);
        xl += dx;
    }

    return Y[0];

}

void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double Y1[256], Y2[256], Y3[256];
    double k1[256], k2[256], k3[256], k4[256];

    // Calcolo di k1
    RHSFunc(t, Y, k1);
    for (k=0;k<neq;k++) Y1[k] = Y[k] + 0.5*dt*k1[k];

    // Calcolo di k2
    RHSFunc(t + 0.5*dt, Y1, k2);
    for (k=0;k<neq;k++) Y2[k] = Y[k] + 0.5*dt*k2[k];

    // Calcolo di k3
    RHSFunc(t + 0.5*dt, Y2, k3);
    for (k=0;k<neq;k++) Y3[k] = Y[k] + dt*k3[k];

    // Calcolo di k4
    RHSFunc(t + dt, Y3, k4);
    for (k=0;k<neq;k++) Y[k] += (dt/6.0) * (k1[k] + 2*k2[k] + 2*k3[k] + k4[k]);

}


void TridiagonalSolver(int n, double a[], double b[], double c[], double r[], double x[]){

    int nmax = 2048;
    if (n>nmax) throw runtime_error("n is too large");
    double h[nmax], p[nmax];
    
    h[0] = c[0]/b[0];
    p[0] = r[0]/b[0];

    for(int i=0;i<n-1;i++){
        h[i+1] = c[i+1]/(b[i+1] - a[i+1]*h[i]);
        p[i+1] = (r[i+1] - a[i+1]*p[i])/(b[i+1] - a[i+1]*h[i]);
    }

    x[n-1] = p[n-1];
    for(int i=n-2;i>=0;i--){
        x[i] = p[i] - h[i]*x[i+1];
    }

}


double Bisection(double (*func)(double), double a, double b, double tol, int& ntry){
    
    double xm, f_xm;         
    double sup_lim, inf_lim;  
    double f_inf, f_sup;   
    const double epsilon = 1.e-12;  // Tolleranza per confronti con zero   

    sup_lim = b;
    inf_lim = a;
    xm = (a + b) / 2.0;
    f_xm = func(xm);
    f_inf = func(a);
    f_sup = func(b);
    ntry = 1;  

    // Verifica se uno degli estremi è già una radice
    if(fabs(f_sup) < epsilon) return sup_lim;
    else if(fabs(f_inf) < epsilon) return inf_lim;

    while(fabs(sup_lim - inf_lim) > tol){

        if(fabs(f_xm) < epsilon) break;  // Trovata una radice esatta

        // Verifica se la radice è nel sottointervallo [inf_lim, xm]
        else if(f_inf * f_xm < 0){
            sup_lim = xm;                
            xm = (sup_lim + inf_lim) / 2.0;  
            f_xm = func(xm);              
        }
        // Altrimenti la radice è nel sottointervallo [xm, sup_lim]
        else if(f_inf * f_xm > 0){
            inf_lim = xm;                 
            f_inf = f_xm;                 
            xm = (sup_lim + inf_lim) / 2.0;  
            f_xm = func(xm);              
        }

        ntry++;  
    }
    return xm;
}