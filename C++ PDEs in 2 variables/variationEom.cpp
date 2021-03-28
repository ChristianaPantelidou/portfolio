#include "variationEom.hpp"
#include <time.h>
#define PI 3.14159265358979
#include <fstream>
#include "auxiliaryFunctions.h"



//function that computes the bc at the UV
void varBcUV(double* result,int i,int j,int i1, int j1, cFunction  * f,double * xgrid, double *ygrid, derivativeClass & diffOperator, int Nf,int Nx,int Ny,double* c){

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
    
   // #include "varBcUVXLine.cpp"
   #include "varBcUV.cpp"
    
}

//function that computes the bc at the U
void varBcIR(double* result,int i,int j,int i1, int j1, cFunction * f, double * xgrid, double *ygrid,derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){

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

    
   // #include "varBcIRXLine.cpp"
     #include "varBcIR.cpp"
    
}


//function that computes the bc at the U
void varBcIRmod(double* result,int i,int j,int i1, int j1, cFunction * f, double * xgrid, double *ygrid,derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){
    
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
    
    
    // #include "varBcIRXLine.cpp"
#include "varBcIRmod.cpp"
    
}


//function that computes the bc at the UV
void varEom(double* result,int i,int j,int i1, int j1, cFunction  * f,double * xgrid, double *ygrid,derivativeClass & diffOperator, int Nf,int Nx,int Ny,double* c){
    
    clock_t t;
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

   // t=clock();
    #include "varEom.cpp"
    //t = clock() - t;
  // std::cout<<((double)t)/CLOCKS_PER_SEC<<"inside";
}


void variationEom(Eigen::SparseMatrix<double> &mat, cFunction *f,double * xgrid, double *ygrid, derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2));
    int i1,j1,i,j,l,f1;
    clock_t t;

    t=clock();
    double * result= new double[Nf*Nf]();
    
    //std::cout<<"here"<<"\n";
    //here i want to rewrite the loops of the variation matrix in the oposite way i did before... this is just to improve redability of the code, no other gain
    //Here i'm focusing on one at the time, and I perturb all the rest
    
        for (int i=1;i<Nx-1; i++){
            for (int j=0;j<Ny-1; j++){//Looping over all points within the domain
                //t=clock();
            
                //Perturb all the points along  a vetical line
                for (int i1=0;i1<Nx; i1++){//eomVar localised vline
                    j1=j;
                    if (i1!=i){
                        varEom(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                        for (int l=0;l<Nf; l++){
                            for (int f1=0;f1<Nf; f1++){
                                 tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                             }
                        }
                    }
                }
                

               //Perturb all the points along  a horizontal line
                for (int j1=0;j1<Ny-1; j1++){//eomVar localised hline
                    i1=i;
                    if (j1!=j){
                        varEom(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                         for (int l=0;l<Nf; l++){
                             for (int f1=0;f1<Nf; f1++){
                                 tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                             }
                        }
                    }
                    
                }
                
                
                //Perturb at the same point
                j1=j;
                i1=i;
                varEom(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                 for (int l=0;l<Nf; l++){
                     for (int f1=0;f1<Nf; f1++){
                         tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                         //std::cout<<result[f1];
                     }
                }

               
            }
        }
    
  //  std::cout<<"here2"<<"\n";
    
   // t=clock();
        for (int j=0;j<Ny-1; j++){//looping over all points in the UV
            i=0;
            
            //Perturb all the points along  a vetical line
            for (int i1=1;i1<Nx; i1++){
                    j1=j;
                    varBcUV(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                    for (int l=0;l<Nf; l++){
                        for (int f1=0;f1<Nf; f1++){
                            //std::cout<<result[f1];
                            tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                           
                        }
                    }
            }
            
            //Perturb all the points along  a horizontal line
            for (int j1=0;j1<Ny-1; j1++){
                i1=i;
                if (j!=j1){
                    varBcUV(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                     for (int l=0;l<Nf; l++){
                         for (int f1=0;f1<Nf; f1++){
                             tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                         }
                    }
                }
            }
            
            
            //Perturb at the same Point
            i1=i;
            j1=j;
            varBcUV(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
             for (int l=0;l<Nf; l++){
                 for (int f1=0;f1<Nf; f1++){
                     //std::cout<<result[f1];
                     tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                 }
            }
            //std::cout<<"here";
        }
   //std::cout<<"here3"<<"\n";
   // t = clock() - t;
   // std::cout<<((double)t)/CLOCKS_PER_SEC<<"loop2";
   //  t=clock();
   
    for (j=1;j<Ny-1; j++){//looping over all points in the IR except the first one
        i=Nx-1;
        //Perturb all the points along  a vetical line
        for (i1=0; i1<Nx-1; i1++){
            j1=j;
            varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
            for (int l=0;l<Nf; l++){
                for (int f1=0;f1<Nf; f1++){
                    tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                }
            }
        }
        
        //Perturb all the points along  a horizontal line
        for (j1=0;j1<Ny-1; j1++){
            i1=i;
            if (j!=j1){
                varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                for (int l=0;l<Nf; l++){
                    for (int f1=0;f1<Nf; f1++){
                        tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                    }
                }
            }
        }
        
        
        //Perturb at the same Point
        i1=i;
        j1=j;
        varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
        for (int l=0;l<Nf; l++){
            for (int f1=0;f1<Nf; f1++){
                tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
            }
        }
        //std::cout<<"here";
    }
    
        j=0;//the first point in the IR
        i=Nx-1;
        //Perturb all the points along  a vetical line
        for (i1=0; i1<Nx-1; i1++){
            j1=j;
            varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
            for (int l=0;l<Nf; l++){
                for (int f1=0;f1<Nf; f1++){
                    tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                }
            }
        }
        
        //Perturb all the points along  a horizontal line
        for (j1=0;j1<Ny-1; j1++){
            i1=i;
            if (j!=j1){
                varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                for (int l=0;l<Nf; l++){
                    for (int f1=0;f1<Nf; f1++){
                        tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                    }
                }
            }
        }
        
        
        //Perturb at the same Point
        i1=i;
        j1=j;
        varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
        for (int l=0;l<Nf; l++){
            for (int f1=0;f1<Nf; f1++){
                tripletList.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
            }
        }
        //std::cout<<"here";
    }
    
    
    //t = clock() - t;
    //std::cout<<((double)t)/CLOCKS_PER_SEC<<"loop3";
    //std::cout<<"here4"<<"\n";
    
    std::cout<<((double)t)/CLOCKS_PER_SEC<<"loop2";
    t=clock();
    
    //Eigen::SparseMatrix<double> mat(Nf*Nx*(Ny-1),Nf*Nx*(Ny-1));//my sparse matrix;
    mat.setFromTriplets(tripletList.begin(), tripletList.end());//construct the matrix from the triplets
    //std::cout<<"here5"<<"\n";
    double * varmat= new double[Nf*Nx*(Ny-1)*Nf*Nx*(Ny-1)]();
    
    for (int i=0;i<Nf*Nx*(Ny-1); i++){
        for (int j=0;j<Nf*Nx*(Ny-1); j++){
            varmat[i*Nf*Nx*(Ny-1)+j]=mat.coeffRef(i,j);
        }
    }
    
    writeInFile("varMattest.bin",varmat,Nf*Nx*(Ny-1)*Nf*Nx*(Ny-1));
   // std::cout<<"here6"<<"\n";
    delete[] result;
   // std::cout<<"here";
    
    
}
