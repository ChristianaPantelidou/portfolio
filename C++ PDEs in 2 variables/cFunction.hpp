#ifndef cFunction_hpp
#define cFunction_hpp

#include "function.h"
#include <string>
#include "DerivativeClass.hpp"
#include "auxiliaryFunctions.h"
#include <stdlib.h>
#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

class cFunction{
private:
    int i,Nx,Ny;
    double* Dx;
    double* Dy;
    functionClass element[6]; // elements
    std::map<std::string,int> stringToInt;
public:
    void create(int,int);
    void set(double*,derivativeClass &);
    void get(double*,int);
    void getd10(double*,int);
    void getd20(double*,int);
    double diff(std::string, int,int);
};

#endif