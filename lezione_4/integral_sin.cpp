#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

double Simpson(double (*F)(double), double a, double b, int N);
double GaussQuadratureRule (double (*F)(double), double a, double b, int N, int Ng);
double Rectangular(double (*F)(double), double a, double b, int N);
double Trapezoidal(double (*F)(double), double a, double b, int N);
double func(double x);


int main(){

    cout << setiosflags(ios::scientific) << setprecision(8);

    cout <<  "Gauss, N=1, Ng=3:   \t" << GaussQuadratureRule(func, 0., 0.8, 1, 3) << endl;
    cout <<  "Trapezoidal, n=1:   \t" << Trapezoidal(func, 0, 0.8, 1 ) <<endl;
    cout <<  "Trapezoidal, n=2:   \t" << Trapezoidal(func, 0, 0.8, 2 ) <<endl;
    cout <<  "Trapezoidal, n=4:   \t" << Trapezoidal(func, 0, 0.8, 4 ) <<endl;
    cout <<  "Trapezoidal, n=8:   \t" << Trapezoidal(func, 0, 0.8, 8 ) <<endl;
    cout <<  "Simpson, n=2:       \t" << Simpson(func, 0., 0.8, 2) << endl;
    cout <<  "Simpson, n=4:       \t" << Simpson(func, 0., 0.8, 4) << endl;
    cout <<  "Simpson, n=8:       \t" << Simpson(func, 0., 0.8, 8) << endl;


    ofstream fdata;         
    fdata.open("integral_sin.dat");

    double integral = 0.;
    double x = 0.;

    while(x < 25.){

        integral += GaussQuadratureRule(func, x, (x+0.1), 1, 3);
        fdata << x << "  " << integral << " " << endl;   
        x += 0.1;

    }

    fdata.close(); 

    return 0;
}

double func(double x){

    if(abs(x) > 1.e-6) return sin(x)/x;
    else return 1. - (x*x)/6.;

 }



double Simpson(double (*F)(double), double a, double b, int N) {

    if (N % 2 != 0) {
        // N deve essere pari per la regola di Simpson
        N++;
    }

    double h = (b-a)/N;
    double sum = F(a) + F(b);  
     
    for (int i=1;i<N;i+=2) {
        sum += 4. * F(a+ i*h);  // Punti dispari
    }
    
    for (int i=2;i<N-1;i+=2) {
        sum += 2. * F(a + i*h);  // Punti pari
    }

    return (h/3.)*sum;
}


double Rectangular (double (*F)(double), double a, double b, int N){  
 
    double h = (b-a)/N;
    double sum = 0.;

    for(int i=0;i<N;i++){
        sum += F(a + i*h);
    }

    return sum*h;
 }

double Trapezoidal (double (*F)(double), double a, double b, int N){  
 
    double h = (b-a)/N;
    double sum = (F(a) + F(b))/2.;

    for(int i=1;i<N;i++){
        sum += F(a + i*h);
    }

    return sum*h;
 }


double GaussQuadratureRule(double (*F)(double), double a, double b, int N, int Ng){

    double w[8];
    double x[8];

    if (Ng == 2) {
        x[0] =  sqrt(1.0 / 3.0);     w[0] = 1;
        x[1] = -sqrt(1.0 / 3.0);     w[1] = 1;
    } 
    else if (Ng == 3) {
        x[0] =  0;						 w[0] = 8.0 / 9.0;
        x[1] =  sqrt(3.0 / 5.0);		 w[1] = 5.0 / 9.0;
        x[2] = -sqrt(3.0 / 5.0);		 w[2] = 5.0 / 9.0;
    }
    else if (Ng == 4){
        x[0] =  sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[0] = (18 + sqrt(30)) / 36;
        x[1] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[1] = (18 + sqrt(30)) / 36;
        x[2] =  sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[2] = (18 - sqrt(30)) / 36;
        x[3] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));	w[3] = (18 - sqrt(30)) / 36;
    }
    else if (Ng == 5){
        x[0] =  0;												w[0] = 128.0 / 225.0;
		x[1] =  1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		w[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0; 
		x[2] = -1.0 / 3.0 * sqrt(5 - 2 * sqrt(10.0 / 7.0));		w[2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
		x[3] =  1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		w[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
		x[4] = -1.0 / 3.0 * sqrt(5 + 2 * sqrt(10.0 / 7.0));		w[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
    }

    double sum = 0.;
    double xc;
    double h = (b-a)/N;
    double sumk;

    for(int i=0;i<N;i++){
        xc = ((a + i*h) + (a + (i+1)*h))/2.;
        sumk = 0.;
        for(int k=0;k<Ng;k++) {
            sumk += (w[k] * F(h*0.5*x[k] + xc));
        } 
        sum += sumk;
    }

    return sum*h*0.5;
    
}

