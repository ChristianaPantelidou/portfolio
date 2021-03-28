#ifndef function_hpp
#define function_hpp

#include <stdlib.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "DerivativeClass.hpp"

class functionClass{
private:
    int Nx,Ny;
    double* functionProfile;
    int m,l;
public:
    void create(int, int);
    double value(int, int);
    void set( double*);
    void setProfile( double*);
    void getProfile(double *,int);
    //double* derivative(std::string, double*, double*);
    double* derivative(std::string,derivativeClass &);
};

#endif