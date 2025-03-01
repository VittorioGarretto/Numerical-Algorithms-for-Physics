#include <iostream> 
using namespace std;

int main() {

    cout << "Give me 2 integers:" << endl;
    int i,j;
    cin >> i >> j;
    cout << endl;

    cout << "Give me 2 floats:" << endl;
    float x,y;
    cin >> x >> y;
    cout << endl;


    cout << "Operations between the integers:" << endl;
    cout << i << " + " << j << " = " << i+j << endl;
    cout << i << " - " << j << " = " << i-j << endl;
    cout << i << " x " << j << " = " << i*j << endl;
    if(j != 0) {

        cout << i << " / " << j << " = " << (float)i/j << endl;

    }
    else {

        cout << "Cannot divide by zero" << endl;

    }


    cout << "Operations between the floats:" << endl;
    cout << x << " + " << y << " = " << x+y << endl;
    cout << x << " - " << y << " = " << x-y << endl;
    cout << x << " x " << y << " = " << x*y << endl;
    if(y != 0) {

        cout << x << " / " << y << " = " << x/y << endl;

    }
    else {

        cout << "Cannot divide by zero" << endl;

    }

    return 0; 
    
}