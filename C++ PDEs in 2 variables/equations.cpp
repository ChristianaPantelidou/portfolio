#include "equations.hpp"
#define PI 3.141592653589793

//function that computes the bc at the UV
void bcUV(double* result, int Nf,int i,int j,cFunction * f,double* xgrid, double* ygrid,double* c){
    double &s=c[0];
    double &n=c[1];
    double &zh=c[2];
    double &beta=c[3];
    double &L=c[4];
    double &cV=c[5];
    double &ct=c[6];
    double &cu=c[7];
    double &cv1=c[8];
    double &cv2=c[9];
    
    cFunction &Qtt=f[0];
    cFunction &Qzz=f[1];
    cFunction &Qxx=f[2];
    cFunction &Qyy=f[3];
    cFunction &Qz1=f[4];
    cFunction &ay=f[5];
    cFunction &by=f[6];
    cFunction &h=f[7];
    
    #include "bcUV.cpp"
    //result[0]=Qtt.diff("d01",i,j);
}


//function that computes the bc at the UV
void bcIR(double* result,int Nf,int i,int j,cFunction * f,double* xgrid, double* ygrid,double* c){
    double &s=c[0];
    double &n=c[1];
    double &zh=c[2];
    double &beta=c[3];
    double &L=c[4];
    double &cV=c[5];
    double &ct=c[6];
    double &cu=c[7];
    double &cv1=c[8];
    double &cv2=c[9];
    
    cFunction &Qtt=f[0];
    cFunction &Qzz=f[1];
    cFunction &Qxx=f[2];
    cFunction &Qyy=f[3];
    cFunction &Qz1=f[4];
    cFunction &ay=f[5];
    cFunction &by=f[6];
    cFunction &h=f[7];
    
    #include "bcIR.cpp"
    //std::cout<<result[0]<<"\n";
   
   // for (int l = 0 ; l<5; l++){
    //    std::cout<<c[l]<<"\n";
    //}
}

//function that computes the bc at the UV
void bcIRmod(double* result,int Nf,int i,int j,cFunction * f,double* xgrid, double* ygrid,double* c){
    double &s=c[0];
    double &n=c[1];
    double &zh=c[2];
    double &beta=c[3];
    double &L=c[4];
    double &cV=c[5];
    double &ct=c[6];
    double &cu=c[7];
    double &cv1=c[8];
    double &cv2=c[9];
    
    cFunction &Qtt=f[0];
    cFunction &Qzz=f[1];
    cFunction &Qxx=f[2];
    cFunction &Qyy=f[3];
    cFunction &Qz1=f[4];
    cFunction &ay=f[5];
    cFunction &by=f[6];
    cFunction &h=f[7];
    
#include "bcIR.cpp"
    //std::cout<<result[0]<<"\n";
    
    // for (int l = 0 ; l<5; l++){
    //    std::cout<<c[l]<<"\n";
    //}
}

//function that computes the equations of motion
void eom(double* result,int Nf,int i,int j,cFunction * f,double* xgrid, double* ygrid,double* c){
    //std::vector result;
    double &s=c[0];
    double &n=c[1];
    double &zh=c[2];
    double &beta=c[3];
    double &L=c[4];
    double &cV=c[5];
    double &ct=c[6];
    double &cu=c[7];
    double &cv1=c[8];
    double &cv2=c[9];
    
    cFunction &Qtt=f[0];
    cFunction &Qzz=f[1];
    cFunction &Qxx=f[2];
    cFunction &Qyy=f[3];
    cFunction &Qz1=f[4];
    cFunction &ay=f[5];
    cFunction &by=f[6];
    cFunction &h=f[7];

    //double test;
    
    //std::cout<<"hereEom";
    
    #include "eom.cpp"
    //std::cout<<Qtt.diff("d00",i,j);
    

}


void totalEqns(double * b ,int Nf,int Nx,int Ny, cFunction* f,double* xgrid, double* ygrid,double* c){//
    int i,j;
    clock_t t;
    
    //t = clock();
    double tt= omp_get_wtime();
    //Parallelising here takes longer: 0.04sec against 0.03sec
#pragma omp parallel
{   double* result= new double[Nf]();
    
    #pragma omp for private(j)
    for (int i = 0 ; i<Nx; i++)
        {
            for (int j = 0 ; j<(Ny-1); j++)
            {
            
                if (i==0){
                    bcUV(result,Nf,i,j,f,xgrid,ygrid,c);//evaluate the bcUV
                   // std::cout<<i<<j<<result[0]<<"\n";
                }else if (i==Nx-1){
                    if (j==0){
                        bcIRmod(result,Nf,i,j,f,xgrid,ygrid,c);
                    }else {
                        bcIR(result,Nf,i,j,f,xgrid,ygrid,c);
                    }
                    //std::cout<<i<<j<<result[0]<<"\n";
                }else {
                    eom(result,Nf,i,j,f,xgrid,ygrid,c);//or eom
                }
                for (int l = 0 ; l<Nf; l++){//and i put them in the matrix M
                    b[l*Nx*(Ny-1)+i*(Ny-1)+j]=result[l];
                //std::cout<<i<<j<<result[l]<<"\n";
                }
            
            }
        }
    delete []result;//this is not needed any more
}
   
    /* double* result= new double[Nf]();
     for (int i = 0 ; i<Nx; i++){
        for (int j = 0 ; j<(Ny-1); j++){// at every cycle, pwEqns hold Nf values
            if (i==0){
                bcUV(result,Nf,i,j,f,xgrid,ygrid,c);//evaluate the bcUV at point i,j
            }else if (i==Nx-1){
                bcIR(result,Nf,i,j,f,xgrid,ygrid,c);//or bcIR
            }else {
                eom(result,Nf,i,j,f,xgrid,ygrid,c);//or eom
            }
            for (int l = 0 ; l<Nf; l++){
                b[l*Nx*(Ny-1)+i*(Ny-1)+j]=result[l];
            }
    
        }
    }
     delete []result;//this is not needed any more*/
    
   // t = clock() - t;
    tt=omp_get_wtime( )-tt;
    std::cout<<((double)tt)<<"_time_eqns"<<"\n";

}
    
