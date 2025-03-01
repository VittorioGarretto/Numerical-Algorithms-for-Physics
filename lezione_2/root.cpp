#include <iostream>
using namespace std;
#include <cmath>
#include <iomanip>

int main(){

    double num, guess, root;

    cout << "Enter real number:" ;
    cin >> num;

    cout << "Enter a guess for the square root:";
    cin >> guess;

    double x =  0.5*(guess + num/guess);

    while(abs(x-guess)>1.e-12){

        guess = x;
        x = 0.5*(guess + num/guess);
        
    }

    root = guess;

    cout << setiosflags(ios::scientific);

    cout << "The sqrt of " << num  << setprecision(14) << " is " << sqrt(num) << setprecision(14) << endl;
    cout << "The iteration sqrt of " << num << setprecision(14) << " is " << root << setprecision(14) << endl;


    return 0;


}