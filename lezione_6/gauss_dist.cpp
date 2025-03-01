#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

double Gauss(double x);

int main(){

    ofstream fdata;         
    fdata.open("gauss_dist.dat");

    double x,y;
    int N = 2500;
    double c = 1.;
    srand48(time(NULL));

    while(N>1){
        x = -5 + 10*drand48(); // numero random in [-5,5]
        y = c*drand48();

        if(y < Gauss(x)){
            fdata << x << " " << y << " " << endl;
            N -= 1;
        }
    }

    fdata.close();

    return 0;
}

double Gauss(double x){

    return (1./(0.5 * sqrt(2*M_PI))) * exp((-1./2.) * (x/0.5) * (x/0.5));

}