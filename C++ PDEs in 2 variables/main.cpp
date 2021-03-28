//boundary conditions

//hard coded spectral, no finite diff

//***to fix*** Limitation:1 function, one equation
//***Hard coded*** number of functions and other parameters appearing in the equation

#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Core>
#include <Accelerate/Accelerate.h>
#include "cFunction.hpp"
#include "DerivativeClass.hpp"
#include "variationEom.hpp"
#include "equations.hpp"
#include "newtonRaphson.hpp"
#include "xiSquare.hpp"
#include<Eigen/SparseLU>
#include <time.h>
#include <limits>
//#include <omp.h>

typedef std::numeric_limits< double > dbl;

using namespace std;

int main(void){
    double domainx=1.;
    double domainy=1.;
    double threshold=0.0000001;
    double epsilon=0.0001;
    double maxStep,maxXi;
    double* c=new double[10]();//A fixed length array to store all the constants that appear in the problem
    double* dim=new double[3]();
    int fdim,ydim,xdim,i,j,m,l,k;
    
   readFromFile("dimensions.bin",dim,3);
    fdim=(int)dim[0];
    xdim=(int)dim[1];
    ydim=(int)dim[2];

    
/////////////////Read file with Paramenters///////////////
    double tt= omp_get_wtime();
    readFromFile("parameters.bin",c,10);
    writeInFile("parameterstext.bin",c,10);
////////////////////Creating the lattice///////////////////////////
 
    //Fourier lines: [0,1]
    //Chebyshev lines: [-1,1]
    double* xgrid = new double[xdim];
    double* ygrid = new double[ydim];
    
   	//Chebyshev grid line
    readFromFile("xgrid.bin",xgrid,xdim);
    writeInFile("xgridtest.bin",xgrid,xdim);
    
    //Fourier grid line
    readFromFile("ygrid.bin",ygrid,ydim);
    writeInFile("ygridtest.bin",ygrid,ydim);

    
////////Create a vector of FunctionClass members, set dimensions and initialise all of them////////////

    derivativeClass diffOperator(xdim,ydim);
    cFunction* f=new cFunction[fdim];//create a dynamical array of cFunctions. Note that the language does not allow a constructor with variables in this case.
   for (int l=0; l<fdim; l++){
       f[l].cFunction::create(xdim,ydim);
   }
    

    double* guessTot = new double[fdim*xdim*(ydim-1)];
    double* guessFun= new double[xdim*(ydim-1)];
    readFromFile("initialGuess.bin",guessTot,fdim*xdim*(ydim-1));//read  guess of each function at a time and sent it over
    writeInFile("initialGuesstest.bin",guessTot,fdim*xdim*(ydim-1));
    
    for (int l=0; l<fdim; l++){
            for (int m=0; m<xdim*(ydim-1); m++){
                guessFun[m]=guessTot[l*xdim*(ydim-1)+m];
            }
            f[l].cFunction::set(guessFun, diffOperator);
    }
    
    delete[] guessTot;
    delete[] guessFun;
    
    double* step = new double[fdim*xdim*(ydim-1)]();
    double* xi = new double[(xdim-2)*(ydim-1)]();
    double* stepFun= new double[xdim*(ydim-1)]();
//Newton-Raphson
    do{
        newtonRaphson(step,f,xgrid,ygrid,diffOperator,fdim,xdim,ydim,c);
        
        
       // t = clock();
        //update cFunction
       
        for (int l=0; l<fdim; l++){
            for (int m=0; m<xdim*(ydim-1); m++){
                stepFun[m]=-step[l*xdim*(ydim-1)+m];
            }
            f[l].cFunction::set(stepFun, diffOperator);
        }
        //std::cout<<"updates";
       // t = clock() - t;
       // std::cout<<((double)t)/CLOCKS_PER_SEC;
 
        computeXi(xi,f,xgrid,ygrid,xdim,ydim,c);
        maxXi=findMax(xdim*(ydim-1),xi);
        
        //calculate maxStep
        maxStep=findMax(fdim*xdim*(ydim-1),step);
        std::cout<<maxStep<<"_MaxStep"<<"\n";
    }while (maxStep>threshold) ;

    //delete[] stepFun;
    //delete[] step;

    
    double* result=new double[fdim*xdim*(ydim-1)]();
    
  for (int i=0;i<fdim;i++){
       f[i].get(result,i);
    }
    writeInFile("solution.bin", result, fdim*xdim*(ydim-1));

    for (int i=0;i<fdim;i++){
        f[i].getd10(result,i);
    }
    writeInFile("solution10.bin", result, fdim*xdim*(ydim-1));

    
    for (int i=0;i<fdim;i++){
        f[i].getd20(result,i);
    }
    writeInFile("solution20.bin", result, fdim*xdim*(ydim-1));
    
    
    computeXi(xi,f,xgrid,ygrid,xdim,ydim,c);
     writeInFile("xiVal.bin", xi, (xdim-2)*(ydim-1));
    
    
    std::cout<<"end";
    delete[] result;
    delete[] stepFun;
    //delete[] step;
    tt=omp_get_wtime( )-tt;
    std::cout<<((double)tt)<<"_time_Total"<<"\n";
 return 0;
}
