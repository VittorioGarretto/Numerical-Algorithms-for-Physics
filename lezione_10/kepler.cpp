#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);

int main(){

    double Y[4];
    double alpha = 0.3;

    Y[0] = 4.;
    Y[1] = 0.;
    Y[2] = 0.;
    Y[3] = sqrt(alpha/4.);

    int neq = 4;
    double dt = 0.1 * sqrt(Y[0]*Y[0] + Y[1]*Y[1])/sqrt(Y[3]*Y[3] + Y[2]*Y[2]);  // 0.1 * r/v
    //double dt = 2*M_PI*20*sqrt(4*4*4)/10000.;
    int counter = 0;
    double vel = Y[2];
    
    ofstream fdata;         
    fdata.open("kepler.dat");
    double t = 0.;
    fdata << Y[0] << " " << Y[1] << endl;
    while(counter < 21){
        RK4Step(t, Y, dYdt, dt, neq);
        t += dt;
        fdata << Y[0] << " " << Y[1] << endl;
        if(vel*Y[2] < 0) counter++;
        vel = Y[2];
        dt = 0.1 * sqrt(Y[0]*Y[0] + Y[1]*Y[1])/sqrt(Y[3]*Y[3] + Y[2]*Y[2]);
    }
    fdata.close();


    ofstream fdata1;         
    fdata1.open("kepler_max.dat");
    t = 0.;
    counter = 0;
    Y[0] = 4.;
    Y[1] = 0.;
    Y[2] = 0.;
    Y[3] = sqrt(alpha/4.);
    vel = sqrt(alpha/4.);
    dt = 0.1 * sqrt(Y[0]*Y[0] + Y[1]*Y[1])/sqrt(Y[3]*Y[3] + Y[2]*Y[2]);
    fdata1 << t << " " << Y[0] << endl;
    while(counter < 21){
        RK4Step(t, Y, dYdt, dt, neq);
        t += dt;
        fdata1 << t << " " << Y[0] << endl;
        if(vel*Y[2] < 0) counter++;
        vel = Y[2];
        dt = 0.1 * sqrt(Y[0]*Y[0] + Y[1]*Y[1])/sqrt(Y[3]*Y[3] + Y[2]*Y[2]);
    }
    fdata1.close();

    return 0;
}

void dYdt(double t, double *Y, double *R){

    double x = Y[0];
    double y = Y[1];
    double vx = Y[2];
    double vy = Y[3]; 

    R[0] = vx;
    R[1] = vy;
    R[2] = -x/((x*x + y*y)*sqrt(x*x + y*y));
    R[3] = -y/((x*x + y*y)*sqrt(x*x + y*y));

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