
#include "xiSquare.hpp"



//function that computes the equations of motion
double pwComputeXi(int i,int j,cFunction * f,double* xgrid, double* ygrid,double* c){
    double result;
    double &s=c[0];
    double &n=c[1];
    double &zh=c[2];
    double &beta=c[3];
    double &L=c[4];
    
    cFunction &Qtt=f[0];
    cFunction &Qzz=f[1];
    cFunction &Qxx=f[2];
    cFunction &Qyy=f[3];
    cFunction &Qz1=f[4];
    cFunction &ay=f[5];
    cFunction &by=f[6];
    cFunction &h=f[7];
    
    
    #include "xi2.cpp"
    
    return result;
    
}

void computeXi(double * b, cFunction* f,double* xgrid, double* ygrid,int Nx, int Ny,double* c){
    int i,j;

    for (int i = 1 ; i<(Nx-1); i++){
        for (int j = 0 ; j<(Ny-1); j++){
            b[i*(Ny-1)+j]=pwComputeXi(i,j,f,xgrid,ygrid,c);
        }
    }

}