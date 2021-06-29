//
//  matrix.hpp
//  
//
//  Created by Christiana Pantelidou on 29/06/2021.
//

#ifndef matrix_hpp
#define matrix_hpp
#include <iostream>
#include <stdio.h>

#endif /* matrix_hpp */

class matrix
{
 private:
  int m;
  int n;
  double **M;


  public:
    matrix();
    matrix(const int& nx, const int& ny);
    matrix(double ** A,const int& nx,const  int& ny);
    matrix(const matrix & A);
    ~matrix();
    double**extractArray();
    int extractM();
    int extractN();
    double& operator() (const int & i, const int& j);
    matrix& operator=(const matrix & A);
    matrix operator+(const matrix& A)const;
    matrix operator+(const double a)const;
    matrix operator-(const matrix& A)const;
    matrix operator-(const double a)const;
    matrix operator*(const matrix & A)const;
    matrix operator*(const double a)const;
    matrix operator/(const double a)const;
    friend std::ostream& operator<<(std::ostream &out,const matrix &A);
};


