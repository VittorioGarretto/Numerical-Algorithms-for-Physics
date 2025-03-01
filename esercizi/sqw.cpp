// s = 1 ;  Root: 7.12337528e+00 ;  ntry: 34
// s = 2 ;  Root: 2.81527537e+01 ;  ntry: 34
// s = 3 ;  Root: 6.14337274e+01 ;  ntry: 34
// ================================
// s = 1 ;  Root: 7.12337527e+00 ;  ntry: 6
// s = 2 ;  Root: 2.81527537e+01 ;  ntry: 5
// s = 3 ;  Root: 6.14337274e+01 ;  ntry: 5
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
using namespace std;

double Secant(double (*func)(double), double a, double b, double tol, int& ntry);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
double func(double x);

int g_s = 1;
double g_V1 = 250.;
double g_V2 = 80.;

int main(){

    double tol = 1.e-8;
    int ntry = 0;
    cout << setiosflags(ios::scientific) << setprecision(8);

    double a = 0.; //dominio per sqrt
    double b = g_V2; // dominio per asin

    for(g_s=1;g_s<=3;g_s++){
        cout << "s = " << g_s << " ; " << " Root: " <<  Bisection(func, a, b, tol, ntry)
        << " ; " << " ntry: " << ntry << endl;
    }

    cout << "================================" << endl;

    g_s = 1;
    a = 0.;
    b = 10.;
    cout << "s = " << g_s << " ; " << " Root: " <<  Secant(func, a, b, tol, ntry)
    << " ; " << " ntry: " << ntry << endl;

    g_s = 2;
    a = 20.;
    b = 30.;
    cout << "s = " << g_s << " ; " << " Root: " <<  Secant(func, a, b, tol, ntry)
    << " ; " << " ntry: " << ntry << endl;

    g_s = 3;
    a = 55.;
    b = 70.;
    cout << "s = " << g_s << " ; " << " Root: " <<  Secant(func, a, b, tol, ntry)
    << " ; " << " ntry: " << ntry << endl;

    return 0;
}

double func(double x){
    
    double s, V1, V2;
    s = g_s;
    V1 = g_V1;
    V2 = g_V2;
    return sqrt(x) - M_PI*s + asin(sqrt(x/V1)) + asin(sqrt(x/V2));

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

double Secant(double (*func)(double), double a, double b, double tol, int& ntry) {

    double sup_lim, inf_lim;   
    double f_inf, f_sup;       
    double dx;                 

    inf_lim = a;
    sup_lim = b;
    dx = sup_lim - inf_lim;
    f_inf = func(inf_lim);
    f_sup = func(sup_lim);
    ntry = 0;                 

    while (fabs(dx) > tol) {
        // Calcolo dell'incremento dx usando la formula del metodo della secante
        dx = f_sup * (sup_lim - inf_lim) / (f_sup - f_inf);
        inf_lim = sup_lim;        
        f_inf = f_sup;            
        sup_lim = sup_lim - dx;   
        f_sup = func(sup_lim);    

        ntry++;                   
    }

    return sup_lim;               
}