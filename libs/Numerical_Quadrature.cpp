#include "Numerical_Quadrature.h"


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

double Simpson(double (*F)(double), double a, double b, int N) {

    if (N % 2 != 0) {
		cout << "Simpson(): Invalid argument: N must be even." << endl;  // N deve essere pari
		exit(1);
	}
	
	const double h = (b-a)/N;		
	double sum = F(a) + F(b);			
	double xi = a;
	double w = 4.0;						

	for (int i=1;i<N;i++) {
		xi += h;
		sum += F(xi) * w;
		w = 6.0 - w;					// Alternanza tra w = 4 e w = 2
	}
	
	return sum*h/3.;
}

double Gauss(double (*F)(double), double a, double b, int N, int Ng){

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
    double xc;  // centro di ciascun sottointervallo
    double h = (b-a)/N;
    double sumk;

    for(int i=0;i<N;i++){  // Ciclo su ciascuno dei N sottointervalli
        xc = ((a + i*h) + (a + (i+1)*h))/2.;
        sumk = 0.;
        for(int k=0;k<Ng;k++) {
            sumk += (w[k] * F(h*0.5*x[k] + xc)); // Regola di gauss sui sottointervalli
        } 
        sum += sumk;
    }

    return sum*h*0.5;
    
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