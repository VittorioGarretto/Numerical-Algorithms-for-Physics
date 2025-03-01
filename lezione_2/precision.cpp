#include <iostream>
using namespace std;
#include <cmath>

int main(){

    float x = 1.;
    float eps = 1.;

    while(true){

        if((x+eps)==x)
            break;
        else
            eps = eps/10.;

    }

    cout << "The precison is: " << eps << endl;

    return 0;

}
