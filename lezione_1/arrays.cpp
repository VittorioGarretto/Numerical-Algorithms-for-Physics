#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

void stat(int, vector<float>&, float&, float&, float&);

int main(){

    int lenght;
    cout << "Define the lenght of the array:" << endl;
    cin >> lenght;
    float mean = 0, std = 0, var = 0;

    vector<float> array(lenght);

    srand48(time(NULL));

    for(int i=0;i<lenght;i++){

        array[i] = drand48();

    }

    stat(lenght, array, mean, var, std);

    cout << endl;
    cout << "The mean of array is: " << mean << endl;
    cout << "The std of array is: " << std << endl;
    cout << "The var of array is: " << var << endl;

    return 0;
}


void stat(int lenght, vector<float>& array, float& mean, float& var, float& std){

    float total = 0;
    
    for(int i=0;i<lenght;i++){

        total += array[i];

    }

    mean = total/lenght;
    total = 0;

    for(int i=0;i<lenght;i++){

        total += pow((array[i]-mean),2);

    }

    var = total/lenght;
    std = pow(var,0.5);

    cout << "The array is:" << endl;
    for(int i=0;i<lenght;i++){

        cout << array[i] << " ";

    }

}

