#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void TridiagonalSolver(int n, double a[], double b[], double c[], double r[], double x[]);

int main(){

    const int n = 5;
    double a[n], b[n], c[n], x[n], r[n];

    a[0] = NAN; a[1] = 1; a[2] = 1; a[3] = 1; a[4] = 1;
    b[0] = 2;   b[1] = 2; b[2] = 2; b[3] = 2; b[4] = 2;
    c[0] = 1;   c[1] = 1; c[2] = 1; c[3] = 1; c[4] = NAN;
    r[0] = 1;   r[1] = 0; r[2] = 3; r[3] = 1; r[4] = 0;

    TridiagonalSolver(n, a, b, c, r, x);

    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

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