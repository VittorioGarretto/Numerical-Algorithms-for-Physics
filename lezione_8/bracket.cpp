#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <stdexcept>
using namespace std;

double func(double x);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
void Bracketing(double (*func)(double), double a, double b, int N, double XL[], double XR[], int& k);

int main(){

    double a = -10.;
    double b = 10.;
    int N = 10;
    double XL[N];
    double XR[N];
    int k = 0;
    int ntry = 0;

    Bracketing(func, a, b, N, XL, XR, k);

    cout << setiosflags(ios::scientific) << setprecision(8);

    for(int j=0;j<k;j++){
        cout << "root " << j+1 << " :" << Bisection(func, XL[j], XR[j], 1.e-8, ntry) << endl;
    }

    return 0;
}

double func(double x){
    return sin(x) - ((x/10.)*(x/10.) + (x/5.) + 1./3.);
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

void Bracketing(double (*func)(double), double a, double b, int N, double XL[], double XR[], int& k) {

    double xl, xr;
    double f_curr, f_prev;
    double h = (b-a)/N;
    k = 0;  // Inizializza il contatore degli intervalli

    f_prev = func(a); // Calcola f(a) una volta

    for (int i = 0; i < N; i++) {
        xl = a + i*h;           
        xr = a + (i + 1) * h;     
        f_curr = func(xr); // Calcola f(xr)

        if (f_prev * f_curr < 0) {
            XL[k] = xl;          
            XR[k] = xr;          
            k++;                            
        }

        f_prev = f_curr; 
    }
}