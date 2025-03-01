#include <iostream> 
using namespace std;

int main() {

    cout << "Print numbers from 1 to 10 using for loop:" << endl;

    for(int i=0;i<10;i++){

        cout << i+1 << " ";

    }

    cout << endl;
    cout << "Print numbers from 1 to 10 using while loop:" << endl;
    int j = 0;

    while(j<10){

        cout << j+1 << " ";
        j++;

    }

    cout << endl;
    cout << "Print odd numbers from 1 to 10:" << endl;

    for(int i=0;i<11;i++){

        if(i%2 != 0)
            cout << i << " ";

    }

    cout << endl;
    cout << "Print even numbers from 1 to 10:" << endl;

    for(int i=0;i<11;i++){

        if(i%2 == 0)
            cout << i << " ";

    }

    cout << endl;
    return 0; 
    
}