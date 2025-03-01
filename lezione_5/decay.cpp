#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

int main(){

    ofstream fdata;         
    fdata.open("decay.dat");

    double lambda = 0.01;
    int N = 1000;
    double prob;
    int N_new;
    int timestep = 0;

    N_new = N;

    srand48(time(NULL));
    while (N > 10){
        for(int i=0;i<N;i++){
            prob = drand48();
            if (prob <= lambda){
                N_new -= 1;
            }
        }

        N = N_new;
        timestep += 1;
        fdata << timestep << " " << (double)N/1000 << " " << endl;
    }

    fdata.close();

    return 0;
}