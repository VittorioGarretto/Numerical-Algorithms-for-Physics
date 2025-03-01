#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double func(double x);
double Rectangular(double (*F)(double), double a, double b, int N);
double Trapezoidal(double (*F)(double), double a, double b, int N);
double Simpson(double (*F)(double), double a, double b, int N);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(8);

    cout <<  "Rectangular: \t" << Rectangular(func, 0, 1, 4 ) <<  endl;
    cout <<  "Trapezoidal: \t" << Trapezoidal(func, 0, 1, 4 ) << endl;
    cout <<  "Simpson:     \t" << Simpson(func, 0, 1, 4) << endl;

    cout << "=================================" << endl;

    double tol = 1.e-5;
    int N_1 = 4;
    

    double I_N = Rectangular(func, 0, 1, N_1);
    double I_half_N = Rectangular(func, 0, 1, N_1/2);
    while(abs(I_N - I_half_N)>tol){

        I_half_N = I_N;
        N_1 *= 2;
        I_N = Rectangular(func, 0, 1, N_1);

    }

    int N_2 = 4;
    I_N = Trapezoidal(func, 0, 1, N_2);
    I_half_N = Trapezoidal(func, 0, 1, N_2/2);
    while(abs(I_N - I_half_N)>tol){

        I_half_N = I_N;
        N_2 *= 2;
        I_N = Trapezoidal(func, 0, 1, N_2);

    }


    int N_3 = 4;
    I_N = Simpson(func, 0, 1, N_3);
    I_half_N = Simpson(func, 0, 1, N_3/2);
    while(abs(I_N - I_half_N)>tol){

        I_half_N = I_N;
        N_3 *= 2;
        I_N = Simpson(func, 0, 1, N_3);

    }

    double sum = 0.;
    cout <<  "Rectangular: \t" << Rectangular(func, 0, 1, N_1) << " N=" << N_1 << endl;
    cout <<  "Trapezoidal: \t" << Trapezoidal(func, 0, 1, N_2) << " N=" << N_2 <<endl;
    cout <<  "Simpson: \t" << Simpson(func, 0, 1, N_3) << " N=" << N_3 <<endl;

    return 0;
}



double func(double x){

    return exp(-x);

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

 /*double Simpson(double (*F)(double), double a, double b, int N){

    double h = (b-a)/N;
    double sum = 0;
    for(int i=0;i<N+1;i++){

        if(i==0 || i==N){
            sum += F(a + i*h) * (h/3.);
        }
        else if(i%2==0){
            sum += F(a + i*h) * (2.*h/3.);
        }
        else if(i%2!=0){
            sum += F(a + i*h) * (4.*h/3.);
        }
    }

    return sum;

 }*/

 
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
 