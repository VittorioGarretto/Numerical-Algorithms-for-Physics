#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>
using namespace std;

void DotProduct(int **A, int *b, int *c, int n);

int main(){

    int nrow = 4;
    int ncol = 4;
    int **A; 
    A    = new int*[nrow]; 
    A[0] = new int [ncol*nrow]; 
    for (int i = 1; i < nrow; i++) A[i] = A[i-1] + ncol;

    A[0][0] = 1;  A[0][1] = 3;  A[0][2] = 2;  A[0][3] = -4;
    A[1][0] = 7;  A[1][1] = 2;  A[1][2] = 4;  A[1][3] = 1;
    A[2][0] = 0;  A[2][1] = -1; A[2][2] = 2;  A[2][3] = 2;
    A[3][0] = 6;  A[3][1] = 3;  A[3][2] = 0;  A[3][3] = 1;

    // Print the matrix  
    for (int i = 0; i < nrow; i++){    
        for (int j = 0; j < ncol; j++){      
            cout << setw(10) << right << A[i][j] << "  ";    
        }    
        cout << endl;  
    }  

    int *b; 
    b = new int[4];  

    b[0] = 1;
    b[1] = 0;
    b[2] = 3;
    b[3] = 2;

    cout << "============================================================" << endl;
    for (int i = 0; i < 4; i++){    
       cout << setw(10) << right << b[i] << "  "; 
       cout << endl;   
    }    
    cout << "============================================================" << endl;

    int *c = new int[4];

    DotProduct(A, b, c, 4);

    for (int i = 0; i < 4; i++){    
       cout << setw(10) << right << c[i] << "  "; 
       cout << endl;   
    } 

    delete[] b;
    delete[] c;
    delete[] A[0];
    delete[] A;

    return 0;
}


void DotProduct(int **A, int *b, int *c, int n){

    int i,j;

    for(i=0;i<n;i++){
        c[i] = 0;
        for(j=0;j<n;j++) c[i] += A[i][j]*b[j];
    }

}