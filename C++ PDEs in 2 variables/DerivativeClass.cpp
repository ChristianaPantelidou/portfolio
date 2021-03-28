#include "DerivativeClass.hpp"
#include <stdlib.h>

derivativeClass::derivativeClass(int xdim,int ydim){
    
    Nx=xdim;
    Ny=ydim;
    
    D0=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    Dx=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    Dy=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    Dx2=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    Dy2=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    Dxy=new double[Nx*(Ny-1)*Nx*(Ny-1)];
    
    readFromFile("D0.bin",D0,Nx*Nx*(Ny-1)*(Ny-1));
    readFromFile("Dx.bin",Dx,Nx*Nx*(Ny-1)*(Ny-1));
    readFromFile("Dy.bin",Dy,Nx*Nx*(Ny-1)*(Ny-1));
    readFromFile("Dx2.bin",Dx2,Nx*Nx*(Ny-1)*(Ny-1));
    readFromFile("Dy2.bin",Dy2,Nx*Nx*(Ny-1)*(Ny-1));
    readFromFile("Dxy.bin",Dxy,Nx*Nx*(Ny-1)*(Ny-1));
    
    writeInFile("Dxtest.bin",Dx,Nx*Nx*(Ny-1)*(Ny-1));
}

double derivativeClass::matrixValue(std::string der,int i,int j){
    double D;
    
    if (der=="d00"){D=D0[i*(Ny-1)*Nx+j];}
    else if (der=="d10"){D=Dx[i*(Ny-1)*Nx+j];}
    else  if (der=="d01"){D=Dy[i*(Ny-1)*Nx+j];}
    else  if (der=="d20"){D=Dx2[i*(Ny-1)*Nx+j];}
    else  if (der=="d02"){D=Dy2[i*(Ny-1)*Nx+j];}
    else  if (der=="d11"){D=Dxy[i*(Ny-1)*Nx+j];}
    
    return D;
}
