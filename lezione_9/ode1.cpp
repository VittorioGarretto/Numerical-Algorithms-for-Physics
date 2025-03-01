#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);

int main(){

    double Y[1];
    Y[0] = 1;
    int neq = 1;
    double dt = 0.5;
    double t = 0.;

    ofstream fdata;         
    fdata.open("ode1.dat");

    cout << setiosflags(ios::scientific) << setprecision(6);

    cout << t << " " << Y[0] << " " <<  fabs(exp(-t*t/2.) - Y[0]) << endl;
    fdata << t << " " << Y[0] << endl;
    while(t<3.){
        EulerStep(t, Y, dYdt, dt, neq);
        t += dt;
        cout << t << " " << Y[0] << " " <<  fabs(exp(-t*t/2.) - Y[0]) << endl;
        fdata << t << " " << Y[0] << endl;
    }
    fdata.close();

    return 0;
}

void dYdt(double t, double *Y, double *R){

    double y = Y[0]; 
    R[0] = -t*y;

}

void EulerStep(double t,  double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq){

    int k; 
    double rhs[256];

    RHSFunc (t, Y, rhs); 
    for (k = 0; k < neq; k++) Y[k] += dt*rhs[k];
}
