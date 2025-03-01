#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

int main(){

    ofstream fdata;         
    fdata.open("pi.dat");

    int N = 10;
    int N_acc = 0;
    double x, y;
    double I;
    srand48(time(NULL));

    while(fabs((I/M_PI) - 1.) > 1.e-4){

        for(int i=0;i<N;i++){
            x = -1 + 2*drand48();
            y = -1 + 2*drand48();
            if(sqrt(x*x + y*y) <= 1) N_acc += 1;
        }

        I = 4. * double(N_acc)/N;
        fdata << N << " " << abs((I/M_PI) - 1.) << " " << endl;
        N_acc = 0;
        N *= 2;
    }

    fdata.close();

    return 0;
}