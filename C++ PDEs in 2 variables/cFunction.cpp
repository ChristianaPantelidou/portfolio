#include "cFunction.hpp"

void cFunction::create(int xdim,int ydim){
    
    Nx=xdim;
    Ny=ydim;
    
   for (int l=0;l<6;l++){
        element[l].functionClass::create(xdim,ydim);
   }
    
    stringToInt["d00"]=0;
    stringToInt["d10"]=1;
    stringToInt["d01"]=2;
    stringToInt["d20"]=3;
    stringToInt["d02"]=4;
    stringToInt["d11"]=5;
    
}

void cFunction::set(double* profile,derivativeClass & diffOperator){
    
    element[0].functionClass::setProfile(profile);
    element[1].functionClass::set(element[0].functionClass::derivative("d10",diffOperator));
    element[2].functionClass::set(element[0].functionClass::derivative("d01",diffOperator));
    element[3].functionClass::set(element[0].functionClass::derivative("d20",diffOperator));
    element[4].functionClass::set(element[0].functionClass::derivative("d02",diffOperator));
    element[5].functionClass::set(element[0].functionClass::derivative("d11",diffOperator));
    
}

void cFunction::get(double* solution,int l){
    
element[0].functionClass::getProfile(solution,l);
    
}

void cFunction::getd10(double* solution,int l){
    
    element[1].functionClass::getProfile(solution,l);
    
}

void cFunction::getd20(double* solution,int l){
    
    element[3].functionClass::getProfile(solution,l);
    
}

double cFunction::diff(std::string der, int i,int j){
    int n;
    double result;
    
    
    n=stringToInt[der];
    result=element[n].functionClass::value(i,j);
    
    return result;
}
