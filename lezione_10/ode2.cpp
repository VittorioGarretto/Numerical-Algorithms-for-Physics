#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);

int main(){

    double Y_eu[2];
    double Y_k2[2];
    double Y_k4[2];

    Y_eu[0] = 1.;
    Y_eu[1] = 0.;
    Y_k2[0] = 1.;
    Y_k2[1] = 0.;
    Y_k4[0] = 1.;
    Y_k4[1] = 0.;

    int neq = 2;
    int nsteps = 200;
    double t = 20. * M_PI;
    double dt = t/double(nsteps);
    
    ofstream fdata;         
    fdata.open("ode2_eu.dat");
    t = 0.;
    fdata << Y_eu[0] << " " << Y_eu[1] << endl;
    for(int i=0;i<nsteps;i++){
        EulerStep(t, Y_eu, dYdt, dt, neq);
        t += dt;
        fdata << Y_eu[0] << " " << Y_eu[1] << endl;
    }
    fdata.close();

    ofstream fdata1;         
    fdata1.open("ode2_k2.dat");
    t = 0.;
    fdata1 << Y_k2[0] << " " << Y_k2[1] << endl;
    for(int i=0;i<nsteps;i++){
        RK2Step(t, Y_k2, dYdt, dt, neq);
        t += dt;
        fdata1 << Y_k2[0] << " " << Y_k2[1] << endl;
    }
    fdata1.close();

    ofstream fdata2;         
    fdata2.open("ode2_k4.dat");
    t = 0.;
    fdata2 << Y_k4[0] << " " << Y_k4[1] << endl;
    for(int i=0;i<nsteps;i++){
        RK4Step(t, Y_k4, dYdt, dt, neq);
        t += dt;
        fdata2 << Y_k4[0] << " " << Y_k4[1] << endl;
    }
    fdata2.close();

    cout << setiosflags(ios::scientific) << setprecision(6);
    cout << t << " " << Y_eu[0] << " " <<  Y_eu[1] << endl;
    cout << t << " " << Y_k2[0] << " " <<  Y_k2[1] << endl;
    cout << t << " " << Y_k4[0] << " " <<  Y_k4[1] << endl;

    //////////////////////////////////////////////////////////////////
    t = 3.;
    nsteps = 2;
    dt = t/double(nsteps);

    ofstream fdata3;         
    fdata3.open("ode2_eu_err.dat");
    while(nsteps < 2050){
        t = 0.;
        Y_eu[0] = 1.;
        Y_eu[1] = 0.;
        for(int i=0;i<nsteps;i++){
            EulerStep(t, Y_eu, dYdt, dt, neq);
            t += dt;
        }
        fdata3 << dt << " " << fabs(Y_eu[0] - cos(t)) << endl;
        nsteps *= 2;
        t = 3.;
        dt = t/double(nsteps);
    }  
    fdata3.close();


    t = 3.;
    nsteps = 2;
    dt = t/double(nsteps);

    ofstream fdata4;         
    fdata4.open("ode2_k2_err.dat");
    while(nsteps < 2050){
        t = 0.;
        Y_k2[0] = 1.;
        Y_k2[1] = 0.;
        for(int i=0;i<nsteps;i++){
            RK2Step(t, Y_k2, dYdt, dt, neq);
            t += dt;
        }
        fdata4 << dt << " " << fabs(Y_k2[0] - cos(t)) << endl;
        nsteps *= 2;
        t = 3.;
        dt = t/double(nsteps);
    }  
    fdata4.close();

    t = 3.;
    nsteps = 2;
    dt = t/double(nsteps);

    ofstream fdata5;         
    fdata5.open("ode2_k4_err.dat");
    while(nsteps < 2050){
        t = 0.;
        Y_k4[0] = 1.;
        Y_k4[1] = 0.;
        for(int i=0;i<nsteps;i++){
            RK4Step(t, Y_k4, dYdt, dt, neq);
            t += dt;
        }
        fdata5 << dt << " " << fabs(Y_k4[0] - cos(t)) << endl;
        nsteps *= 2;
        t = 3.;
        dt = t/double(nsteps);
    }  
    fdata5.close();


    return 0;
}

void dYdt(double t, double *Y, double *R){

    double x = Y[0];
    double y = Y[1]; 
    R[0] = y;
    R[1] = -x;

}

void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double rhs[256];

    RHSFunc (t, Y, rhs); 
    for (k = 0; k < neq; k++) Y[k] += dt*rhs[k];
}

void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double rhs[256];
    double Y1[256];
    double k1[256];
    double k2[256];

    RHSFunc(t, Y, k1);
    for (k = 0; k < neq; k++) {
        Y1[k] = Y[k] + 0.5 * dt * k1[k];
    }

    // Calcolo di k2
    RHSFunc(t + 0.5 * dt, Y1, k2);

    // Aggiornamento finale di Y
    for (k = 0; k < neq; k++) {
        Y[k] += dt * k2[k];
    }

}

void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double rhs[256];
    double Y1[256], Y2[256], Y3[256];
    double k1[256], k2[256], k3[256], k4[256];

    // Calcolo di k1
    RHSFunc(t, Y, k1);
    for (k = 0; k < neq; k++) {
        Y1[k] = Y[k] + 0.5 * dt * k1[k];
    }

    // Calcolo di k2
    RHSFunc(t + 0.5 * dt, Y1, k2);
    for (k = 0; k < neq; k++) {
        Y2[k] = Y[k] + 0.5 * dt * k2[k];
    }

    // Calcolo di k3
    RHSFunc(t + 0.5 * dt, Y2, k3);
    for (k = 0; k < neq; k++) {
        Y3[k] = Y[k] + dt * k3[k];
    }

    // Calcolo di k4
    RHSFunc(t + dt, Y3, k4);

    // Aggiornamento finale di Y
    for (k = 0; k < neq; k++) {
        Y[k] += (dt / 6.0) * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
    }

}

