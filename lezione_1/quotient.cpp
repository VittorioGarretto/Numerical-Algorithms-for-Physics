#include <iostream>
using namespace std;

void Quotient(int, int);

int main(){

    int a,b;
    cout << "Give me 2 integers:" << endl;
    cin >> a >> b;

    Quotient(a, b);

    return 0;
}

void Quotient(int a, int b){

    cout << a << " / " << b << " = " << a/b << endl;
    cout << "the remainder is " << a%b << endl;

}