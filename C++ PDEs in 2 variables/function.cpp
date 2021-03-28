#include "function.h"


//Set dimension, create the functionProfile with the right dimesnion and reads the matrices
//This could not have been done with a constructor because parameters should be sent over and this can not be dome when you call the class with dynamical allocation
void functionClass::create(int xdim, int ydim){
    
    Nx=xdim;
    Ny=ydim;
    functionProfile=new double[Nx*(Ny-1)]();
  }


//updates the value of the profile to be something
void functionClass::setProfile(double* profile){
    
  for (int m=0; m<Nx*(Ny-1); m++){
      //std::cout<<"here";
      functionProfile[m]=functionProfile[m]+profile[m];
    }
 
}
//sets the value of the profile to be something, to initialise
void functionClass::set(double* profile){
    
    for (int m=0; m<Nx*(Ny-1); m++){
        //std::cout<<"here";
        functionProfile[m]=profile[m];
    }
    
}

//returns the value of the function at a given point
double functionClass::value(int i,int j){
    return functionProfile[i*(Ny-1)+j];
}

//returns the value of the function at a given point
void functionClass::getProfile(double *result,int l){
    for (int i=0;i<Nx*(Ny-1);i++){
         result[l*Nx*(Ny-1)+i]=functionProfile[i];
    }
}
//Computes the derivative of the function at all points and returns a vector derivativeProfile
double* functionClass::derivative(std::string der, derivativeClass &diffOperator){
    double* derivativeProfile=new double[Nx*(Ny-1)] ();
    
    for (int i=0; i<Nx*(Ny-1);i++){
        for (int l=0; l<Nx*(Ny-1);l++){
            derivativeProfile[i]+=diffOperator.matrixValue(der,i,l)*functionProfile[l];
        }
    }
    return derivativeProfile;
}


