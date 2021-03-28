#include <stdlib.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

class gridLine{
	private:
		std::string gridType;
		int numberOfDimensions;
		Eigen::VectorXd line;
	public:
		gridLine(int,std::string);
		double pwGetValue(int i);
        Eigen::VectorXd getLine();
		std::string getGridType();
		int getNumberOfDimensions();
};