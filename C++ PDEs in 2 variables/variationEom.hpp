#ifndef variationEom_hpp
#define variationEom_hpp

#include <Eigen/Sparse>
#include "cFunction.hpp"
#include "DerivativeClass.hpp"
#include <Eigen/SparseCore>
#include <stdlib.h>
# include "umfpack.h"
//#include <omp.h>

//void variationEom(Eigen::SparseMatrix<double>&, cFunction* ,double *, double *,derivativeClass &,int,int,int,double*);

void variationEom(long *,long *,double *, cFunction* ,double *, double *,derivativeClass &,int,int,int,double*);

#endif