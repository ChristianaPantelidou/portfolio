//
//  matrix.cpp
//  
//
//  Created by Christiana Pantelidou on 29/06/2021.
//

#include "matrix.hpp"
#include <stdio.h>


matrix::matrix(){//default constructor
    m=1;
    n=1;
    M=0;
}
    
matrix::matrix(const int& nx, const int& ny){//constructor
    m=nx;
    n=ny;
    M=new  double *[m];
    for (int i=0;i<m;i++){
        M[i]=new double[n](); //initialised to zero
    }
            
}
    
matrix::matrix(double ** A,const int& nx,const  int& ny){//constructor overload
    m=nx;
    n=ny;
    
    M=new  double *[m];
    for (int i=0;i<m;i++){
        M[i]=new double[n];
    }
    
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            M[i][j]=A[i][j];
        }
    }
}
    
matrix::matrix(const matrix & A){//copy constructor
    m=A.m;
    n=A.n;
    
    M=new  double *[m];
    for (int i=0;i<m;i++){
        M[i]=new double[n];
    }
    
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            M[i][j]=A.M[i][j];
        }
    }
}
        
double** matrix::extractArray(){
    return M;
}
    
int matrix::extractM(){
    return m;
}
int matrix::extractN(){
    return n;
}
    
//        double extractelem(const int & i,const int & j){//extract (i,j) element
//            return M[i][j];
//        }
    
double& matrix::operator() (const int & i, const int& j){//extract (i,j) element
    return M[i][j];
}
    
//
//        void update(const int &i,const int& j,const  double a){//update the (i,j) element
//            M[i][j]=a;
//        }
        
    
//        void update( double ** a){//update the whole matrix, assumes right dimensions
//            M=a;
//        }
//
matrix& matrix::operator=(const matrix & A){//assignement operator
    if (this==&A)
        return *this;
    
    for (int i=0;i<m;i++){
        delete [] M[i];
    }
    M=0;
//
    n=A.n;
    m=A.m;
    M=new  double *[m];
    for (int i=0;i<m;i++){
        M[i]=new double[n];
    }
//
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            M[i][j]=(A.M)[i][j];
        }
    }
        
    return *this;
}
    
matrix matrix::operator+(const matrix& A)const {//add a matrix
    
    
    if (n!=A.n|| m!=A.m){
        throw std::runtime_error("Dimensions do not match");
        return *this;
    }
    else
    {
        matrix C(m,n);
        for (int i=0;i<m;i++){
            for (int j=0;j<n;j++){
                C.M[i][j]=M[i][j]+A.M[i][j];
            }
        }
        return C;
    }

}
//
matrix matrix::operator+(const double a)const{//broadcast the addition of a number
    matrix C(m,n);
            
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            C.M[i][j]=M[i][j]+a;
        }
    }
return C;
}
//
matrix matrix::operator-(const matrix& A)const{//subtract a matrix
            
    if (n!=A.n|| m!=A.m){
        throw std::runtime_error("Dimensions do not match");
        matrix C(m, n);
        return *this;
    }
    else
    {
        matrix C(m,n);
        for (int i=0;i<m;i++){
            for (int j=0;j<n;j++){
                C.M[i][j]=M[i][j]-A.M[i][j];
            }
        }
        return C;
    }
}
    
matrix matrix::operator-(const double a)const{//broadcast the subtraction of a number
    matrix C(m,n);
        
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            C.M[i][j]=M[i][j]-a;
        }
    }
            
    return C;
}
    
    
    
matrix matrix::operator*(const matrix & A)const{//multiplication of 2 matrices
    int n1=A.n;
        
    if (n!=A.m){
        throw std::runtime_error("Dimensions do not match");
        return *this;
    }
    else
    {
        matrix C( m, n1);
        for (int i=0;i<m;i++){
            for (int j=0;j<n1;j++){
                C.M[i][j]=0;
                for (int l=0;l<n;l++){
                    C.M[i][j]+=M[i][l]*(A.M)[l][j];
                }
            }
        }
        return C;
    }
           
}
    
matrix matrix::operator*(const double a)const{//broadcast multiplication with a number
    matrix C(m, n);
            
    for (int i=0;i<m;i++){
        for (int j=0;j<n; j++){
            C.M[i][j]=M[i][j]*a;
        }
    }
    return C;
        
}
    
matrix matrix::operator/(const double a)const{//broadcast division with a number
    matrix C(m, n);
            
    for (int i=0;i<m;i++){
        for (int j=0;j<n; j++){
            C.M[i][j]=M[i][j]/a;
        }
    }
    return C;
        
}
    
 

matrix::~matrix(){//Destructor
    if(M!=0){
        for (int i=0;i<m;i++){
            delete [] M[i];
        }
    }
    M=0;
    m=0;
    n=0;
}
    
std::ostream& operator<<(std::ostream &out,const matrix &A){//print
    out << "{";
    for (int i=0; i<A.m; i++) {
        out<< "{";
        for (int j=0; j<A.n; j++) {
            out<<(A.M)[i][j];
            if (j<A.n-1) {out << ", ";}
            }
        if (i<A.m-1) {out << "},";}
    }
    out<< "}}" << std::endl;
    return out;
}


