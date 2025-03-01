#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
double Residual (double s);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);

int main(){

    double Y[2];
    int neq = 2;
    double dr = 20./1000.;
    double r = 0.;
    double s = 0.2;

    ofstream fdata;         
    fdata.open("poisson_s.dat");
    for (int k=0;k<6;k++){
        Y[0] = 0.;
        Y[1] = k*s;
        fdata << r << " " << Y[0] << endl;
        for(int i=0;i<1000;i++){
            RK4Step(r, Y, dYdt, dr, neq);
            r += dr;
            fdata << r << " " << Y[0] << endl;
        }
        fdata << endl << endl;
        r = 0.;
    }
    fdata.close();

    ////////////////////////////////////////////////////////////

    ofstream fdata1;         
    fdata1.open("poisson_res.dat");
    s = 0.;
    while(s<5.){
        fdata1 << s << " " << Residual(s) << endl;
        s += 0.1;
    }
    
    fdata1.close();

    ///////////////////////////////////////////////////////////

    double s1 = 0.3;
    double s2 = 0.8;
    int ntry = 0;
    double tol = 1.e-6;
    cout << setiosflags(ios::scientific) << setprecision(6);

    double s_true = Bisection(Residual, s1, s2, tol, ntry);
    cout << "Bisection -> s = " << s_true << endl;

    ///////////////////////////////////////////////////////////

    ofstream fdata2;         
    fdata2.open("poisson_phi.dat");
    r = 0.;
    Y[0] = 0.;
    Y[1] = s_true;
    fdata2 << r << " " << Y[0] << endl;
    for(int i=0;i<1000;i++){  // non si devono usare loop su double ma solo su int 
        RK4Step(r, Y, dYdt, dr, neq);
        r += dr;
        fdata2 << r << " " << Y[0] << endl;
    }
    
    return 0;
}

void dYdt(double r, double *Y, double *R){ 

    R[0] = Y[1];
    R[1] = -4. * M_PI * r * (1./(8.*M_PI)) * exp(-r);

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

double Residual(double s){

    double Y[2];
    int neq = 2;
    int nstep = 1000;
    double r = 0.;
    double dr = 20./nstep;
    Y[0] = 0.;
    Y[1] = s;

    for(int i=0;i<nstep;i++){
        RK4Step(r, Y, dYdt, dr, neq);
        r += dr;
    }

    return Y[0] - 1.;

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