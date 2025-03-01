#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void TridiagonalSolver(int n, double a[], double b[], double c[], double r[], double x[]);

int main(){

    const int n = 32;   
    double k = M_PI*M_PI;
    double dx = 0.5/(n-1);   // eq diff y'' = -k y + 1 with y(0) = y(0.5) = 0
    double x0 = 0.;
    double x1 = 0.;
    double a[n], b[n], c[n], x[n], r[n];
    x[0] = x0;
    x[n-1] = x1;

    for(int i=0;i<n;i++) a[i] = 1.;
    a[0] = a[1] = NAN;

    for(int i=0;i<n;i++) b[i] = -2. + k*dx*dx;
    b[0] = b[n-1] = NAN;

    for(int i=0;i<n;i++) c[i] = 1.;
    c[n-1] = c[n-2] = NAN;

    for(int i=0;i<n;i++) r[i] = dx*dx;
    r[1] = dx*dx - x0;
    r[n-2] = dx*dx - x1;
    r[0] = r[n-1] = NAN;

    TridiagonalSolver(n-2, a+1, b+1, c+1, r+1, x+1);

    ofstream fdata;
    fdata.open("bvp_linear.dat");
    for (int i = 0; i < n; i++) {
        fdata << i*dx << " " << x[i] << " " << -1./k * (cos(sqrt(k) * i*dx) + sin(sqrt(k) * i*dx) - 1.) << endl;
    }
    fdata.close();

    return 0;
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



