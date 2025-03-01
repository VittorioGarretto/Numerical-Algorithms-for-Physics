#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

void func1(float, float&, float&, float&);

int main(){

    cout << setiosflags(ios::scientific) << setprecision(6);
    cout << "Compute sqrt(x^2 + 1) - x for large x" << endl;
    cout << string(60, '=') << endl;
    cout << left << setw(15) << "x" 
         << setw(15) << "fx1" 
         << setw(15) << "fx2" 
         << setw(15) << "f(taylor)" 
         << endl;
    cout << string(60, '-') << endl;

    float x1, fx1, fx2, f_taylor;
    x1 = 1.e4;

    for(int i=0;i<7;i++){

        func1(x1, fx1, fx2, f_taylor);
        x1 *= 10.;
    }

    return 0;

}

void func1(float x, float& fx1, float& fx2, float& f_taylor){

    fx1 = sqrt(x*x + 1) - x;
    fx2 = 1./(sqrt(1.+(x*x)) + x);
    f_taylor = 1./(2.*x);

    cout << left << setw(15) << x
         << setw(15) << fx1 
         << setw(15) << fx2 
         << setw(15) << f_taylor
         << endl;

}
