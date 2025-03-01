#include "Ode_Solver.h"

void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double rhs[256];

    RHSFunc (t, Y, rhs); 
    for (k=0;k<neq;k++) Y[k] += dt*rhs[k];

}

void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double Y1[256];
    double k1[256];
    double k2[256];

    RHSFunc(t, Y, k1);
    for (k=0;k<neq;k++) Y1[k] = Y[k] + 0.5*dt *k1[k];

    // Calcolo di k2
    RHSFunc(t + 0.5*dt, Y1, k2);
    for (k=0;k<neq;k++) Y[k] += dt*k2[k];

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

void pVerlet(double *Y, double *v, void (*Acceleration)(double *, double *), double dt, int neq) {

    int n;
    double a[256]; // Array per accelerazioni

    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

    Acceleration(Y, a);

    for(n=0;n<neq;n++) v[n] += dt*a[n];
    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

}