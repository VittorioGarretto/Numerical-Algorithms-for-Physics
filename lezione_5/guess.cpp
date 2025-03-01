#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

int main(){

    srand(time(NULL)); 
    int number = rand()%100+1;

    int guess;
    int iter = 0;
    int i = 0;
    int j = 100;

    while (guess != number){

        cout << "n in [" << i  << "," << j << "]" << endl;
        cout << "type your guess #" << iter << endl;
        cin >> guess;

        if(guess > number){
            j = guess;
            iter += 1;
        }
        else if(guess < number){
            i = guess;
            iter += 1;
        }
        else if(guess == number){
            cout << "Challenge complete" << endl;
        }
    } 
     

    return 0;
}