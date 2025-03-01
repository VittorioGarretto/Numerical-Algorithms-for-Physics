#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void pVerlet(double *Y, double *v, void (*Acceleration)(double *, double *), double dt, int neq);
void Acceleration(double *Y, double *a);
void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);

int main(){

    double Y[2];

    Y[0] = 1.;
    Y[1] = 0.;

    int neq = 2;
    double dt = 0.02 * (2*M_PI);
    double E = 0.5 * Y[1]*Y[1] + 0.5 * Y[0]*Y[0];
    double t = 0.;

    ofstream fdata;         
    fdata.open("harmonic_k2.dat");
    fdata << t << " " << E << endl;
    while(t < 60.){
        RK2Step(t, Y, dYdt, dt, neq);
        t += dt;
        E = 0.5 * Y[1]*Y[1] + 0.5 * Y[0]*Y[0];
        fdata << t << " " << E << endl;
    }
    fdata.close();


    ofstream fdata1;         
    fdata1.open("harmonic_verlet.dat");
    double v[1];
    double Y_pV[1];
    Y_pV[0] = 1.;
    v[0] = 0.;
    neq = 1;
    E = 0.5 * v[0]*v[0] + 0.5 * Y_pV[0]*Y_pV[0];
    t = 0.;
    fdata1 << t << " " << E << endl;
    while(t < 60.){
        pVerlet(Y_pV, v, Acceleration, dt, neq);
        t += dt;
        E = 0.5 * v[0]*v[0] + 0.5 * Y_pV[0]*Y_pV[0];
        fdata1 << t << " " << E << endl;
    }
    fdata1.close();


    return 0;
}

void dYdt(double t, double *Y, double *R){

    double x = Y[0];
    double vx = Y[1]; 
    
    R[0] = vx;
    R[1] = -x;

}

void Acceleration(double *Y, double *a){

    a[0] = -Y[0];

}


void RK2Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double Y1[256];
    double k1[256];
    double k2[256];

    RHSFunc(t, Y, k1);
    for (k=0;k<neq;k++) Y1[k] = Y[k] + 0.5*dt*k1[k];

    // Calcolo di k2
    RHSFunc(t + 0.5*dt, Y1, k2);
    for (k=0;k<neq;k++) Y[k] += dt*k2[k];

}

void pVerlet(double *Y, double *v, void (*Acceleration)(double *, double *), double dt, int neq){

    int n;
    double a[256]; // Array per accelerazioni

    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

    Acceleration(Y, a);

    for(n=0;n<neq;n++) v[n] += dt*a[n];
    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

}



