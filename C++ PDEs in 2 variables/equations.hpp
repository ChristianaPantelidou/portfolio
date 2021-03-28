#ifndef equations_hpp
#define equations_hpp

#include "cFunction.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <cmath>
//#include <omp.h>

void totalEqns(double*,int,int,int, cFunction*,double*, double*,double*);
void eom(double* ,int,int,int,cFunction*,double*, double*,double*);
void bcUV(double* ,int,int,int,cFunction*,double*, double*,double*);
void bcIR(double* ,int,int,int,cFunction*,double*, double*,double*);
void bcIRmod(double* ,int,int,int,cFunction*,double*, double*,double*);

#endif