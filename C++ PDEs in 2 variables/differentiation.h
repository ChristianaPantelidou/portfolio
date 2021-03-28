#include <stdlib.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

class Differentiation{
private:
    int numberOfDimensions;
    std::string gridType;
    std::string gridPeriodic;
    Eigen::MatrixXd matrix;
    Eigen::MatrixXd matrix2;
public:
    Differentiation(int dim,  std::string type, std::string periodic, gridLine * lattice);
    Eigen::MatrixXd GetMatrix();
    Eigen::MatrixXd GetMatrix2();
    double pwGetMatrix(int jx,int jy);
};