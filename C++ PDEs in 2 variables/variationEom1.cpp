#include "variationEom.hpp"
#include <time.h>
#define PI 3.141592653589793
#include <fstream>



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
#include "varBcIR.cpp"
    
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
   // t = clock() - t;
  // std::cout<<((double)t)/CLOCKS_PER_SEC<<"inside";
}


//void variationEom(Eigen::SparseMatrix<double> &mat, cFunction *f,double * xgrid, double *ygrid, derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){

void variationEom(long *Ap,long *Ai,double*Ax, cFunction *f,double * xgrid, double *ygrid, derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){
    
    typedef Eigen::Triplet<double,long> T;
    std::vector<T> tripletList;
    tripletList.reserve(Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2));
    int i1,j1,i,j,l,f1;
    //clock_t t;
    

    //t=clock();
    double tt= omp_get_wtime();
    
#pragma omp parallel
{
    std::vector<T> tripletList_private;
    double * result= new double[Nf*Nf]();
    #pragma omp for private(j,i1,j1,l,f1)
    for ( i=1;i<Nx-1; i++)
        {
            for (j=0;j<Ny-1; j++)
            {//Looping over all points within the domain
                //t=clock();
            
                //Perturb all the points along  a vetical line
                for ( i1=0;i1<Nx; i1++){//eomVar localised vline
                    j1=j;
                    if (i1!=i){
                        //std::cout<<"here";
                        //t=clock();
                        varEom(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                        for (int l=0;l<Nf; l++){
                            for (int f1=0;f1<Nf; f1++){
                                 tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                             }
                        }
                    }
                }
                
               //Perturb all the points along  a horizontal line
                for (j1=0;j1<Ny-1; j1++){//eomVar localised hline
                    i1=i;
                    if (j1!=j){
                        varEom(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                         for (int l=0;l<Nf; l++){
                             for (int f1=0;f1<Nf; f1++){
                                 tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                                 // std::cout<<result[f1];
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
                         tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                         //std::cout<<result[f1];
                     }
                }

            }
        }
    

    
        #pragma omp for private(i,i1,j1,l,f1)
        for ( j=0;j<Ny-1; j++){//looping over all points in the UV
            i=0;
            
            //Perturb all the points along  a vetical line
            for ( i1=1;i1<Nx; i1++){
                    j1=j;
                    varBcUV(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                    for (int l=0;l<Nf; l++){
                        for (int f1=0;f1<Nf; f1++){
                            //std::cout<<result[f1];
                            tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                           
                        }
                    }
            }
            
            //Perturb all the points along  a horizontal line
            for ( j1=0;j1<Ny-1; j1++){
                i1=i;
                if (j!=j1){
                    varBcUV(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                     for (int l=0;l<Nf; l++){
                         for (int f1=0;f1<Nf; f1++){
                             tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
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
                     tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                 }
            }
            //std::cout<<"here";
        }

        #pragma omp for private(i,i1,j1,l,f1)
        for (j=0;j<Ny-1; j++){//looping over all points in the IR
            i=Nx-1;
            //Perturb all the points along  a vetical line
            for (i1=0; i1<Nx-1; i1++){
                j1=j;
                if (j==0){
                    varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                }else {
                    varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                }
               
                 for (int l=0;l<Nf; l++){
                     for (int f1=0;f1<Nf; f1++){
                         tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                     }
                }
            }
            
            //Perturb all the points along  a horizontal line
            for (j1=0;j1<Ny-1; j1++){
                i=Nx-1;
                i1=Nx-1;
                if (j!=j1){
                    if (j==0){
                        varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                    }else {
                        varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
                    }
                     for (int l=0;l<Nf; l++){
                         for (int f1=0;f1<Nf; f1++){
                             tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                         }
                    }
                }
            }
            
            
            //Perturb at the same Point
            i1=Nx-1;
            j1=j;
            if (j==0){
                varBcIRmod(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
            }else {
                varBcIR(result,i,j,i1,j1,f,xgrid,ygrid,diffOperator, Nf,Nx,Ny,c);
            }
             for (int l=0;l<Nf; l++){
                 for (int f1=0;f1<Nf; f1++){
                     tripletList_private.push_back(T(l*Nx*(Ny-1)+i*(Ny-1)+j,f1*Nx*(Ny-1)+i1*(Ny-1)+j1,result[l*Nf+f1]));
                 }
            }
             //std::cout<<"here";
        }
   
    
    // One tread at a time will enter this part of the pragma
    #pragma omp critical
    tripletList.insert(tripletList.end(), tripletList_private.begin(), tripletList_private.end());
    
    delete[] result;
    tripletList_private.erase(tripletList_private.begin(), tripletList_private.end());
}
    
    long *Ti=new long[Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2)];
    long *Tj=new long[Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2)];
    double *Tx=new double[Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2)];
    for (i=0;i<Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2); i++){
        Ti[i]=tripletList[i].row();
        Tj[i]=tripletList[i].col();
        Tx[i]=tripletList[i].value();
    }

    int status;//converts triplets to compressed form, ready for input in UMFpack
    status=umfpack_dl_triplet_to_col(Nf*Nx*(Ny-1), Nf*Nx*(Ny-1), Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2),Ti,Tj,Tx,Ap,Ai,Ax,NULL);
    delete[] Ti,Tj,Tx;

    
   // t = clock() - t;
    tt=omp_get_wtime( )-tt;
    std::cout<<((double)tt)<<"_time_VarMat"<<"\n";

  /*  Eigen::SparseMatrix<double> mat(Nf*Nx*(Ny-1),Nf*Nx*(Ny-1));
    mat.reserve(Eigen::VectorXi::Constant(Nf*Nx*(Ny-1),Nf*(Nx+Ny-2)));//Show many elements per row will be non-zero
    mat.setFromTriplets(tripletList.begin(), tripletList.end());//construct the matrix from the triplets
   // std::cout<<"here";
    
    
    double * varmat= new double[Nf*Nx*(Ny-1)*Nf*Nx*(Ny-1)]();
    
    for (int i=0;i<Nf*Nx*(Ny-1); i++){
        for (int j=0;j<Nf*Nx*(Ny-1); j++){
            varmat[i*Nf*Nx*(Ny-1)+j]=mat.coeffRef(i,j);
        }
    }
    
    writeInFile("varMattest.bin",varmat,Nf*Nx*(Ny-1)*Nf*Nx*(Ny-1));*/
    tripletList.erase(tripletList.begin(), tripletList.end());

    
}
