#ifndef DerivativeClass_hpp
#define DerivativeClass_hpp

#include <string>
#include <stdio.h>
# include "auxiliaryFunctions.h"

class derivativeClass{
private:
    int Nx,Ny;
    double* D0;
    double* Dx;
    double* Dy;
    double* Dx2;
    double* Dy2;
    double* Dxy;
public:
    derivativeClass(int, int);
    //~derivativeClass();
   // double* matrix(std::string);
    double matrixValue(std::string,int,int);
};

#endif
