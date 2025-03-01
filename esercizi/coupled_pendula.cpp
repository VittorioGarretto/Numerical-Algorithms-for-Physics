#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void Acceleration(double *Y, double *a);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
void pVerlet(double *Y, double *v, void (*Acceleration)(double *, double *), double dt, int neq);

int main(){

    double Y_RK4[4], Y_pV[2], v_pV[2];

    int neq = 4;
    Y_RK4[0] = -M_PI/4.;    Y_pV[0] = -M_PI/4.;
    Y_RK4[1] = 0.;          v_pV[0] = 0.;
    Y_RK4[2] = 0.;          Y_pV[1] = 0.;
    Y_RK4[3] = 0.;          v_pV[1] = 0.;
    double dt = 0.1;
    double tol = 0.25;
    double t = 0.;
    double eps1 = fabs(Y_pV[0] - Y_RK4[0])/M_PI;
    double eps2 = fabs(Y_pV[1] - Y_RK4[2])/M_PI;
    int turn_point1 = 0;  
    int turn_point2 = 0;
    int nstep = 0;
    double theta1 = Y_RK4[0];
    double theta2 = Y_RK4[2];

    while(eps1<=tol && eps2<=tol){
        RK4Step(t, Y_RK4, dYdt, dt, neq);
        pVerlet(Y_pV, v_pV, Acceleration, dt, 2);
        t += dt;
        nstep ++;
        eps1 = fabs(Y_pV[0] - Y_RK4[0])/M_PI;
        eps2 = fabs(Y_pV[1] - Y_RK4[2])/M_PI;
        if(theta1 * Y_RK4[0] <= 0.) turn_point1 ++;
        if(theta2 * Y_RK4[2] <= 0.) turn_point2 ++;
        theta1 = Y_RK4[0];
        theta2 = Y_RK4[2];
    }


    cout << "Loop break at nstep = " << nstep << " t = " << t << endl;
    cout << "eps1 = " << eps1 << " eps2 = " << eps2 << endl;
    cout << "tp1 = " << turn_point1 << " tp2 = " << turn_point2 << endl;


    t = 0.;
    Y_RK4[0] = -M_PI/4.;    Y_pV[0] = -M_PI/4.;
    Y_RK4[1] = 0.;          v_pV[0] = 0.;
    Y_RK4[2] = 0.;          Y_pV[1] = 0.;
    Y_RK4[3] = 0.;          v_pV[1] = 0.;
    nstep = 500;

    ofstream fdata;         
    fdata.open("coupled_pendula.dat");
    fdata << t << " " << Y_RK4[0] << " " << Y_RK4[2] << " " << Y_pV[0] << " " << Y_pV[1] << endl;
    for(int i=0;i<nstep;i++){
        RK4Step(t, Y_RK4, dYdt, dt, neq);
        pVerlet(Y_pV, v_pV, Acceleration, dt, 2);
        t += dt;
        fdata << t << " " << Y_RK4[0] << " " << Y_RK4[2] << " " << Y_pV[0] << " " << Y_pV[1] << endl;
    }
    fdata.close();

    return 0;
}

void dYdt(double t, double *Y, double *R){

    double k = 0.8;
    double g = 9.8;

    double theta1 = Y[0];
    double theta2 = Y[2];
    double omega1 = Y[1];
    double omega2 = Y[3];

    R[0] = omega1;
    R[1] = -g*sin(theta1) - k*(sin(theta1) - sin(theta2));
    R[2] = omega2;
    R[3] = -g*sin(theta2) + k*(sin(theta1) - sin(theta2));

}

void Acceleration(double *Y, double *a){

    double k = 0.8;
    double g = 9.8;

    double theta1 = Y[0];
    double theta2 = Y[1];

    a[0] = -g*sin(theta1) - k*(sin(theta1) - sin(theta2));
    a[1] = -g*sin(theta2) + k*(sin(theta1) - sin(theta2));

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
    double a[64]; // Array per accelerazioni

    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

    Acceleration(Y, a);

    for(n=0;n<neq;n++) v[n] += dt*a[n];
    for(n=0;n<neq;n++) Y[n] += 0.5*dt*v[n];

}
