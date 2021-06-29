//
//  matrix.cpp
//
//

#include <iostream>
using namespace std;
#include "matrix.hpp"

int main()
{
    int m1=2;
    int n1=2;
    
    double ** N=new  double *[m1];
    for (int i=0;i<m1;i++){
        N[i]=new double[n1];
    }
    
    for (int i=0;i<m1;i++){
        for (int j=0;j<n1;j++){
            N[i][j]=i;
        }
    }

    matrix A(m1,n1);
    matrix B(N,m1,n1);
    matrix C(A);
    
    matrix m(2, 2);
    m(0,1)=100;

    cout<< m;
    A=(B+B);
    cout<< B;
    cout<< A;

    return 0;
}


