#include <iostream>
using namespace std;

double Sum(double, double);

int main(){

    double a,b;
    cout << "Give me 2 doubles:" << endl;
    cin >> a >> b;

    double c;
    c = Sum(a,b);

    cout << "c = a + b = " << c << endl;

    return 0;
}

double Sum(double a, double b){

    return a+b;

}