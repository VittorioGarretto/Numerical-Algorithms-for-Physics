#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void dYdt(double t, double *Y, double *R);
void RK4Step(double t, double *Y, void (*RHSFunc)(double, double *, double *), double dt, int neq);
double Residual (double k);
double Bisection(double (*func)(double), double a, double b, double tol, int& ntry);
void Bracketing(double (*func)(double), double a, double b, int N, double XL[], double XR[], int& k);

double g_k = 1.;

int main(){

    double Y[2];
    int neq = 2;
    Y[0] = 0.;
    Y[1] = 1.;
    int nstep = 100;
    double dx = 1./100.;
    double x = 0.;

    ofstream fdata;         
    fdata.open("wave_1.dat");
    fdata << x << " " << Y[0] << endl;
    for(int i=0;i<nstep;i++){
        RK4Step(x, Y, dYdt, dx, neq);
        x += dx;
        fdata << x << " " << Y[0] << endl;
    }
    fdata.close();

    /////////////////////////////////////////////

    x = 0.;
    ofstream fdata1;         
    fdata1.open("wave_k.dat");
    for (int j=1;j<6;j++){
        Y[0] = 0.;
        Y[1] = 1.;
        g_k = double(j);
        fdata1 << x << " " << Y[0] << endl;
        for(int i=0;i<nstep;i++){
            RK4Step(x, Y, dYdt, dx, neq);
            x += dx;
            fdata1 << x << " " << Y[0] << endl;
        }
        fdata1 << endl << endl;
        x = 0.;
    }
    fdata1.close();

    ////////////////////////////////////////////

    double k1 = 3.;
    double k2 = 4.;
    int ntry = 0;
    double tol = 1.e-6;
    cout << setiosflags(ios::scientific) << setprecision(6);

    double k_true = Bisection(Residual, k1, k2, tol, ntry);
    cout << "Bisection -> k = " << k_true << " ntry: " << ntry << endl;
    g_k = k_true;

    /////////////////////////////////////////////

    int N = 10;
    double XL[N];
    double XR[N];
    int h = 0;
    ntry = 0;

    Bracketing(Residual, 1., 20., N, XL, XR, h);

    for(int j=0;j<h;j++){
        cout << "root " << j+1 << " :" << Bisection(Residual, XL[j], XR[j], 1.e-8, ntry) << endl;
    }

    ///////////////////////////////////////////////

    Y[0] = 0.;
    Y[1] = 1.;
    x = 0.;
    g_k = 1.884975e+01;

    ofstream fdata2;         
    fdata2.open("wave_ktrue.dat");
    fdata2 << x << " " << Y[0] << endl;
    for(int i=0;i<nstep;i++){
        RK4Step(x, Y, dYdt, dx, neq);
        x += dx;
        fdata2 << x << " " << Y[0] << endl;
    }
    fdata2.close();


    return 0;
}

void dYdt(double r, double *Y, double *R){ 

    R[0] = Y[1];
    R[1] = -Y[0] * g_k*g_k;

}

double Residual(double k){

    double Y[2];
    int neq = 2;
    int nstep = 100;
    double x = 0.;
    double dx = 1./nstep;
    Y[0] = 0.;
    Y[1] = 1.;
    g_k = k;

    for(int i=0;i<nstep;i++){
        RK4Step(x, Y, dYdt, dx, neq);
        x += dx;
    }

    return Y[0];

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

void Bracketing(double (*func)(double), double a, double b, int N, double XL[], double XR[], int& k) {

    double xl, xr;
    double f_curr, f_prev;
    double h = (b-a)/N;
    k = 0;  // Inizializza il contatore degli intervalli

    f_prev = func(a); 

    for (int i = 0; i < N; i++) {
        xl = a + i*h;           
        xr = a + (i + 1)*h;     
        f_curr = func(xr); 

        if (f_prev * f_curr < 0) {
            XL[k] = xl;          
            XR[k] = xr;          
            k++;                            
        }

        f_prev = f_curr; 
    }
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
