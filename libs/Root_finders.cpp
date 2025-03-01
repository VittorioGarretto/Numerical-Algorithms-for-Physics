#include "Root_finders.h"

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


double FalsePos(double (*func)(double), double a, double b, double tol, int& ntry){

    double x1, x2;
    double f_a, f_b;
    double inf_lim, sup_lim;
    double m, c;
    const double epsilon = 1.e-12;  // Tolleranza per confronti con zero 

    inf_lim = a;
    sup_lim = b;
    f_a = func(a);
    f_b = func(b);
    ntry = 1;

    // Controllo delle condizioni iniziali
    if (f_a * f_b > 0) {
        throw invalid_argument("La funzione non cambia segno nell'intervallo fornito.");
    }
    if (fabs(f_b) < epsilon) return sup_lim;
    if (fabs(f_a) < epsilon) return inf_lim;

    x1 = a;
    x2 = b;

    while (fabs(x2 - x1) > tol) {
        f_a = func(inf_lim);
        f_b = func(sup_lim);
        m = (f_b - f_a) / (sup_lim - inf_lim);
        
        // Evitare la divisione per zero
        if (fabs(m) < epsilon) {
            throw runtime_error("Divisione per zero nel calcolo della pendenza.");
        }

        c = f_a - m*inf_lim;
        x1 = x2;
        x2 = -c/m;
        ntry++;

        // Aggiornamento dei limiti in base al segno della funzione
        if (func(x2) * f_a < 0) sup_lim = x2;
        else inf_lim = x2;
    }

    return x2;
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

double Newton(double (*func)(double), double (*dfunc)(double), double a, double b, double tol, int& ntry) {

    double sup_lim = b, inf_lim = a;   
    double xc = (sup_lim + inf_lim) / 2.0;  
    double f = func(xc);
    double f_der = dfunc(xc);

    double x_guess = xc - f/f_der;   // Prima stima della radice
    double dx = x_guess - xc;
    ntry = 1;  

    while (fabs(func(x_guess)) > tol) { 
        xc = x_guess;                  
        f = func(xc);         
        f_der = dfunc(xc); 
        dx = -f / f_der;               
        x_guess = xc + dx;
        ntry ++;
    }

    return x_guess;  
}

double Horner_pol(const double x, const double a[], const int degree, double& dpdx) {
	
    double p = a[degree];
	dpdx = 0.0;
	for (int i = degree-1; i >= 0; i--) {
		dpdx = dpdx * x + p;
		p = p * x + a[i];
	}
	return p;

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