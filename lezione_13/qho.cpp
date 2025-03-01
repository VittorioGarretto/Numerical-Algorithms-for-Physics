#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double x, double *Y, double *R);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
double Residual(double E);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);

double g_E = 0.5;

int main(){

    double Y[2];
    int neq = 2;
    double x = 10.;
    Y[0] = exp(-x*x*0.5);
    Y[1] = -x * Y[0];
    int nstep = 800;
    double dx = -20./double(nstep);

    ofstream fdata;         
    fdata.open("qho_guess.dat");
    fdata << x << " " << Y[0] << endl;
    for(int i=0;i<nstep;i++){
        RK4Step(x, Y, dYdt, dx, neq);
        x += dx;
        fdata << x << " " << Y[0] << endl;
    }
    fdata.close();

    ////////////////////////////////////////

    ofstream fdata1;         
    fdata1.open("qho_res.dat");
    nstep = 1000;
    double dE = 5./double(nstep);
    double E = 0.;
    for(int j=0;j<nstep;j++){
        fdata1 << E << " " << Residual(E) << endl;
        E += dE;
    }
    
    fdata1.close();

    ////////////////////////////////////////

    double E1 = 0.;
    double E2 = 1.;
    int ntry = 0;
    double tol = 1.e-6;
    cout << setiosflags(ios::scientific) << setprecision(6);

    double E_true = Bisection(Residual, E1, E2, tol, ntry);
    cout << "Bisection -> E = " << E_true << " ntry: " << ntry << endl;
    g_E = E_true;

    return 0;
}

void dYdt(double x, double *Y, double *R){ 

    R[0] = Y[1];
    R[1] = -2. * (g_E - 0.5*x*x) * Y[0];

}

double Residual(double E){

    double YR[2], YL[2];
    int neq = 2;
    double xr = -10.; double xl = 10.;
    YR[0] = exp(-xr*xr*0.5); YL[0] = exp(-xl*xl*0.5);
    YR[1] = -xr * YR[0]; YL[1] = -xl * YL[0];
    int nstep = 800;
    double xm = 0.05;
    double dxr = fabs(-10. - xm)/double(nstep);
    double dxl = -(10 - xm)/double(nstep);
    g_E = E;
    
    for(int i=0;i<nstep;i++){
        RK4Step(xr, YR, dYdt, dxr, neq);
        xr += dxr;
    }

    for(int i=0;i<nstep;i++){
        RK4Step(xl, YL, dYdt, dxl, neq);
        xl += dxl;
    }

    return -(YL[1]*YR[0] - YR[1]*YL[0])/sqrt(YL[1]*YR[0]*YL[1]*YR[0] + YR[1]*YL[0]*YR[1]*YL[0]);

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

