#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

void quadratic(double, double);
void quadratic_corrected(double, double);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(6);

    // Intestazione della tabella
    cout << left << setw(15) << "True x1" 
         << setw(15) << "True x2" 
         << setw(15) << "Pred x1" 
         << setw(15) << "Pred x2" 
         << endl;
    cout << string(60, '-') << endl;

    quadratic(2.,-3.);
    quadratic(1.e-5, 1.e8);
    quadratic(1.e-12, 1.e12);


    cout << left << setw(15) << "True x1" 
         << setw(15) << "True x2" 
         << setw(15) << "Pred x1" 
         << setw(15) << "Pred x2" 
         << endl;
    cout << string(60, '-') << endl;

    quadratic_corrected(2.,-3.);
    quadratic_corrected(1.e-5, 1.e8);
    quadratic_corrected(1.e-12, 1.e12);

    return 0;

}

void quadratic(double x1, double x2){

    double a = 1.;
    double b = -(x1 + x2);
    double c = x1 * x2;

    double discriminant = b * b - 4 * a * c;

    double x1_pred = (-b + sqrt(discriminant)) / (2 * a);
    double x2_pred = (-b - sqrt(discriminant)) / (2 * a);

    cout << left << setw(15) << x1 
         << setw(15) << x2 
         << setw(15) << x1_pred 
         << setw(15) << x2_pred 
         << endl;

}

void quadratic_corrected(double x1, double x2){

    double a = 1.;
    double b = -(x1 + x2);
    double c = x1 * x2;
    double x1_pred;
    double x2_pred;

    double discriminant = b * b - 4 * a * c;

    if(b >= 0){
        x1_pred = (-b - sqrt(discriminant)) / (2 * a);
        x2_pred = (2 * c) / (-b - sqrt(discriminant));
    }

    else{
        x1_pred = (2 * c) / (-b + sqrt(discriminant));
        x2_pred = (-b + sqrt(discriminant)) / (2 * a);
    }

    cout << left << setw(15) << x1 
         << setw(15) << x2 
         << setw(15) << x1_pred 
         << setw(15) << x2_pred 
         << endl;


}


