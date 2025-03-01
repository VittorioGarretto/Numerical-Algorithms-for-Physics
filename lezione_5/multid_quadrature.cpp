#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;


double Gaussian2DQuad (double (*Func)(double, double), double xa, double xb, double ya, double yb, int N, int Ng);
double func(double x, double y);
double func2(double x, double y);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(8);

    cout <<  "Gauss_multid: \t" << Gaussian2DQuad(func, -1., 1., -1., 1, 3, 4) << endl;
    cout <<  "Exact:         \t" << 412. / 45. <<endl;

    cout << "=============================" << endl;

    int N = 1;
    double x = Gaussian2DQuad(func2, -1., 1., -1., 1., N, 4);

    while(fabs(x - M_PI) > 1.e-5){

        N += 1;
        x = Gaussian2DQuad(func2, -1., 1., -1., 1, N, 4);

    }
    cout <<  "N:         \t" << N <<endl;



    return 0;
}


double func(double x, double y){

    return (x*x*x*x)*y*y + 2*(y*y)*x*x - y*x*x + 2.;

}

double func2(double x, double y){

    if (sqrt(x*x + y*y) <= 1) return 1.;
    else return 0.;


}



double Gaussian2DQuad (double (*Func)(double, double), double xa, double xb, double ya, double yb, int N, int Ng){

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


    double dx = (xb - xa)/ N;
    double dy = (yb - ya)/ N;

    double sum = 0.0;
    double sumy, sumx;
    double xc, yc;

    for(int j=0;j<N;j++){ // Loop on sub-intervals in y
        for(int i=0;i<N;i++){  // Loop on sub-intervals in x
            
            xc = ((xa + (j+1)*dx) + (xa + j*dx)) / 2.;  // Center of interval (x)
            yc = ((ya + (i+1)*dy) + (ya + i*dy)) / 2.;  // Center of interval (y)
            sumy = 0;

            for(int jk=0;jk<Ng;jk++){
                sumx = 0.0;   // Initialize sum for this interval
                for(int ik=0;ik<Ng;ik++){  // Apply Gaussian rule to sub-interval
                    sumx += w[ik]*Func(xc + 0.5*dx*x[ik], yc + 0.5*dy*x[jk]);
                }
                sumx *= 0.5*dx; 
                sumy += w[jk]*sumx; 
            }

        sumy *= 0.5*dy; 
        sum += sumy;
        }
    }
    
    return sum;
}