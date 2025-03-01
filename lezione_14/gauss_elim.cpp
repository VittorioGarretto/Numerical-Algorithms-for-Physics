#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

int main(){

    const int n = 4;
    double A[n][n] = { {1, 2, 1, -1}, {3, 2, 4, 4}, {4, 4, 3, 4}, {2, 0, 1, 5} };
    double b[n] = {5, 16, 22, 15};
    double g;

    // Eliminazione di Gauss
    for (int k = 0; k < n - 1; k++) { // Loop over the Gk's
        for (int i = k + 1; i < n; i++) {  // Loop over rows
            g = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++)  A[i][j] -= g*A[k][j];
            A[i][k] = 0.0; 
            b[i] -= g * b[k];
        }
    }

    double x[n];

    for (int i = n - 1; i >= 0; i--) {
        double tmp = b[i];
        for (int j = n-1; j > i; j--)
            tmp -= x[j] * A[i][j];
        x[i] = tmp / A[i][i];
    }

    cout << "Soluzioni:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}