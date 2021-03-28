
#ifndef newtonRaphson_hpp
#define newtonRaphson_hpp

#include "cFunction.hpp"
# include "umfpack.h"

# include "auxiliaryFunctions.h"
#include <Eigen/Sparse>
#include "variationEom.hpp"
#include "equations.hpp"
# include "umfpack.h"
#include <Eigen/IterativeLinearSolvers>

void newtonRaphson(double*,cFunction*,double*,double*, derivativeClass &,int,int,int,double*);

#endif