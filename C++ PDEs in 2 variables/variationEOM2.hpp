#ifndef variationEOM2_hpp
#define variationEOM2_hpp

#include <Eigen/Sparse>
#include "cFunction.hpp"
#include "DerivativeClass.hpp"
#include <stdlib.h>
//#include <omp.h>

void variationEom(Eigen::SparseMatrix<double>&, cFunction* ,double *, double *,derivativeClass &,int,int,int,double*);

#endif
