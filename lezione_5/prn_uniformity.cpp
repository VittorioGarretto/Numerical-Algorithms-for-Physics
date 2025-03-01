#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

int main(){

    ofstream fdata;         
    fdata.open("prn_uniformity.dat");

    srand48(time(NULL));
    for (int i=0;i<1000;i++){
        fdata << i << " " << drand48() << " " << endl;
    }

    fdata.close();

    //////////////////////////////////

    ofstream fdata1;
    ofstream fdata2;         
    fdata1.open("prn_uniformity_k1.dat");
    fdata2.open("prn_uniformity_k2.dat");

    double mom1 = 0., mom2 = 0.;
    int N = 10;
    
    while (N<1e6){
        for (int i=0;i<N;i++){
            double x = drand48();
            mom1 += x;;
            mom2 += x*x;
        }
        mom1 /= N;
        mom2 /= N;
        fdata1 << N << " " << fabs(mom1 - (1./2.)) << " " << endl;
        fdata2 << N << " " << fabs(mom2 - (1./3.)) << " " << endl;
        N *= 2;

        mom1 = 0.;
        mom2 = 0.;
    }

    fdata1.close();
    fdata2.close();

    return 0;
}