#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double Simpson(double (*F)(double), double a, double b, int N);
double GaussQuadratureRule (double (*F)(double), double a, double b, int N, int Ng);
double func(double x);
double func2(double x);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(8);

    cout <<  "Simpson: \t" << Simpson(func, 0., 3., 2) << endl;
    cout <<  "Gauss:   \t" << GaussQuadratureRule(func, 0, 3., 1., 3) << endl;
    
    cout << "=======================================" << endl;

    cout <<  "Simpson: \t" << Simpson(func2, -1., 5., 2) << endl;
    cout <<  "Gauss:   \t" << GaussQuadratureRule(func2, -1., 5., 1, 4) << endl;
    cout <<  "Exact:   \t" << -66. / 5. << endl;




    return 0;
}

 

double func(double x){

    return sqrt(1.+x);

 }

 double func2(double x){

    return (1 - x + 2.*x*x +(x*x*x)/2. + (x*x*x*x)/4. - (x*x*x*x*x)/8);

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


