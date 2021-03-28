
//

#ifndef __PDEs__NewtonMethod__
#define __PDEs__NewtonMethod__

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <complex>

#include "function.h"
#include "CustomFunctions.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <petscksp.h>
#include <petscerror.h>
#include <petscdm.h>
#include <petscdmda.h>

typedef std::complex<double> Complex;

struct Triplet {
    MatStencil row;
    MatStencil col;
    PetscScalar val;
};

// CALL: UpdateFunctions<MyType,MyType>(equations, bconditions, BulkLinearOp, BoundaryLinearOp, Q, D, rgrid, xgrid, ygrid,16, ps );
//CALL:  void UpdateFunctions<double, double>( void (*eom),void (*b_cond),void (*l_eom),void (*l_b_cond),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid, grid<double>& Ygrid,My_Int N_functions, double* params){
    

//CALL: ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, Ygrid, N_functions,params, lindices, i_map , tole );
//CALLED THIS WAY: ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, Ygrid, N_functions,params, lindices, i_map , tole );

//Construct the linear operator
template <typename Matrix, typename Vec, typename dType,typename Type>
void ConstructLinearOp(Matrix &TriplVec, Vec &V ,void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k)
					   ,void (*b_cond) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k)
					   ,void (*l_eom)(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
                       ,void (*l_b_cond)(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<dType, Type>* f, derivativeCol<dType>& D, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid ,My_Int N_functions,dType* params, std::vector<My_Int>& Rindices, PetscInt* map, const dType tolerance){

	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();
	My_Int nvars = Ny*Nx+ (2*Rgrid.MaxOrder()+1)*(Nx+Ny)+N_functions*(2*Rgrid.MaxOrder()+1+Nx+Ny);
	My_Int nbdof = Rindices.size()*N_functions;//R=Reduced
	long long res_estimate= (long long)nvars*(long long)nbdof;
    
    int numtasks;
    int rank;
    
//    double omp_start, omp_end;
    
    //earlier in main: MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_size(PETSC_COMM_WORLD, &numtasks);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


#pragma omp parallel
{
	Type *telem1;
	Type *telem3;
	
	telem1 = new Type[N_functions*N_functions];
	telem3 = new Type[N_functions];

	Matrix TempTriplVec;
    //more efficient when you reserve the memory needed
    //to be reserved on each thread
    TempTriplVec.reserve((long long)res_estimate/((long long)omp_get_num_threads()));


	My_Int Rlower;
	My_Int Rupper;
    My_Int Rupper1;
	
	My_Int Ylower;
	My_Int Yupper;
    My_Int Yupper1;
	
    My_Int i1;
    My_Int i2;
    My_Int i3;
    
#pragma omp for
	for (My_Int ii=0; ii<Rindices.size() ; ii++) { //ii controls the point we are looking at

        
        i3= Rindices[ii]/((Np-1)*Nx); // i1, i2 and i3 are labels of the point we are looking at
        i2= (Rindices[ii]%((Np-1)*Nx))/(Np-1);
        i1=1+(Rindices[ii]%((Np-1)*Nx))%(Np-1);//never equal to 0
        //the r indices are flattened first, then the x and then the y
        
      
       /* Rlower=Rgrid.Lower(i1);
        Rupper=Rgrid.Upper(i1)+1;
        
        Rlower= (Rlower>=1)?Rlower:1;
        Rupper= (Rupper<=Np)?Rupper:Np;*/
        
        Ylower=Ygrid.Lower(i3);
        Yupper=Ygrid.Upper(i3);
        
        Ylower= (Ylower>=0)?Ylower:0;
        Yupper= (Yupper<=Ny)?Yupper:Ny;
        
       // std::cout<<f[0].diff(1);
        
        Rgrid.NewtonLowerUpper(i1,Rlower,Rupper);
        
        Rlower= (Rlower>=1)?Rlower:1;
        Rupper= (Rupper<=Np)?Rupper:Np;
        
   	   // Ygrid.NewtonLowerUpper(i3,Ylower,Yupper);
        
       // Ylower= (Ylower>=0)?Ylower:0;
       // Yupper= (Yupper<=Ny)?Yupper:Ny;
       
        if(i3==0){Yupper=Ny;}	//These two lines impose the boundary conditions on the y coordinate
        if(i3==Ny-1){Ylower=0;}	//The variation of the equations set at 0 and Ny-1 depend on the points of the opposite boundaries
        
        
        
        
        
        if(i1<Np-1){ //if  looking in the bulk
            
                    //i=i1
                    for (My_Int j=0; j<Nx; j++) {
                    	if(j==i2){continue;}// continue: Causes the remaining portion of the enclosing for, range-for, while or do-while loop body to be skipped.
                        //So that the lines are not computed more than once!
                        for (My_Int k=Ylower; k<Yupper; k++) {
                            if(k==i3){continue;}
                            l_eom(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i1,j,k);// first point is where we are looking, second where we perturb
                            
                            for (My_Int l3=0; l3<N_functions; l3++) {
                                for(My_Int l=0; l<N_functions; l++){
                                    
                                    if( (modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance)) //&& j!=i2 && k!=i3
                                    {
                                        //#pragma omp critical
                                         TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[k*Nx*(Np-1)*N_functions +j*(Np-1)*N_functions +(i1-1)*N_functions +l3] ,telem1[l*N_functions+l3]));
                                    }
                                    
                                }}}}
                
                //j=i2
                for (My_Int i=Rlower; i< Rupper; i++) {
                    for (My_Int k=Ylower; k<Yupper; k++) {
                            if(k==i3){continue;}
//                            l_eom(telem1, D,f, Rgrid, params,i,i2,k,i1,i2,i3);
                        l_eom(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i,i2,k);
                        
                            for (My_Int l3=0; l3<N_functions; l3++) {
                                for(My_Int l=0; l<N_functions; l++){
                                    
                                    if( modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance) //&& k!=i3
                                    {
                                        //#pragma omp critical
                                        TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[k*Nx*(Np-1)*N_functions +i2*(Np-1)*N_functions +(i-1)*N_functions +l3 ],telem1[l*N_functions+l3]));
                                    }
                                    
                                }}}}
                
                //k=i3
                for (My_Int i=Rlower; i< Rupper; i++) {
                    for (My_Int j=0; j<Nx; j++) {
                        
//                        omp_start = omp_get_wtime();
                        
//                            l_eom(telem1, D,f, Rgrid, params,i,j,i3,i1,i2,i3);
                        l_eom(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i,j,i3);

//                        omp_end = omp_get_wtime();
//                        
//                        std::cout<< i << " "<< i1 << " " << j << " "<< i2 << " "<< omp_end-omp_start <<'\n';
                        
                            for (My_Int l3=0; l3<N_functions; l3++) {
                                for(My_Int l=0; l<N_functions; l++){
                                    
                                    if(modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance )
                                    {
                                        //#pragma omp critical
                                        TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[i3*Nx*(Np-1)*N_functions +j*(Np-1)*N_functions +(i-1)*N_functions +l3 ],telem1[l*N_functions+l3]));
                                    }
                                    
                                }}}}
        }
        
        if (i1==Np-1) { //looking at the boundary, excluding the "same line" cases
         
            //i=i1
            for (My_Int j=0; j<Nx; j++) {
                if(j==i2){continue;}
                for (My_Int k=Ylower; k<Yupper; k++) {
                    if(k==i3){continue;}
                    //                            l_eom(telem1, D,f, Rgrid, params,i1,j,k,i1,i2,i3);
                    l_b_cond(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i1,j,k);
                    
                    for (My_Int l3=0; l3<N_functions; l3++) {
                        for(My_Int l=0; l<N_functions; l++){
                            
                            if( (modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance)) //&& j!=i2 && k!=i3
                            {
                                //#pragma omp critical
                                TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[k*Nx*(Np-1)*N_functions +j*(Np-1)*N_functions +(i1-1)*N_functions +l3 ],telem1[l*N_functions+l3]));
                            }
                            
                        }}}}
            
            //j=i2
            for (My_Int i=Rlower; i< Rupper; i++) {
                for (My_Int k=Ylower; k<Yupper; k++) {
                    if(k==i3){continue;}
                    //                            l_eom(telem1, D,f, Rgrid, params,i,i2,k,i1,i2,i3);
                    l_b_cond(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i,i2,k);
                    
                    for (My_Int l3=0; l3<N_functions; l3++) {
                        for(My_Int l=0; l<N_functions; l++){
                            
                            if( modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance) //&& k!=i3
                            {
                                //#pragma omp critical
                                TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[k*Nx*(Np-1)*N_functions +i2*(Np-1)*N_functions +(i-1)*N_functions +l3 ],telem1[l*N_functions+l3]));

                            }
                            
                        }}}}
            
            //k=i3
            for (My_Int i=Rlower; i< Rupper; i++) {
                for (My_Int j=0; j<Nx; j++) {//x is always pseudospectral
                    
                    //                        omp_start = omp_get_wtime();
                    
                    //                            l_eom(telem1, D,f, Rgrid, params,i,j,i3,i1,i2,i3);
                    l_b_cond(telem1, D,f, Rgrid, Xgrid, Ygrid, params,i1,i2,i3,i,j,i3); //telem1 is NfXNf
                    
                    //                        omp_end = omp_get_wtime();
                    //
                    //                        std::cout<< i << " "<< i1 << " " << j << " "<< i2 << " "<< omp_end-omp_start <<'\n';
                    
                    for (My_Int l3=0; l3<N_functions; l3++) {
                        for(My_Int l=0; l<N_functions; l++){
                            
                            if(modulus<dType, Type>(telem1[l*N_functions+l3])>tolerance )
                            {
                                //#pragma omp critical
                                TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l],map[i3*Nx*(Np-1)*N_functions +j*(Np-1)*N_functions +(i-1)*N_functions +l3 ],telem1[l*N_functions+l3]));
                            }
                            
                        }}}}
            
        }
    }
    
        
    
    
	#pragma omp for
    for (My_Int ii=0; ii<Rindices.size() ; ii++) {
        i3= Rindices[ii]/((Np-1)*Nx);
        i2= (Rindices[ii]%((Np-1)*Nx))/(Np-1);
        i1=1+(Rindices[ii]%((Np-1)*Nx))%(Np-1);// i1 is non-zero
        
        if(i1<Np-1){// if not at the horizon
            eom(telem3,f, Rgrid, Xgrid, Ygrid, params, i1, i2, i3);} //telem3 is a Nf vector
        else if (i1==Np-1) {//if at the horizon
            b_cond(telem3, f, Rgrid, Xgrid, Ygrid, params, i1, i2, i3);}
        else{std::cout << "Ops!" << std::endl;}
		
        for(My_Int l=0; l<N_functions; l++){
				V[map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l]]=telem3[l];
            }
    }
    
	
	delete[] telem1;
	delete[] telem3;
    
	#pragma omp critical
	TriplVec.insert(TriplVec.end(),TempTriplVec.begin(),TempTriplVec.end());
    TempTriplVec.erase(TempTriplVec.begin(),TempTriplVec.end());
	}
}


template<typename dType, typename Type>
void UpdateFunctions( void (*eom) (Type* elem, Cfunction<dType,Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*b_cond) (Type* elem, Cfunction<dType,Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3,dType* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_eom)(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_b_cond)(Type *elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, grid<dType>& grid2, grid<dType>& grid3, dType* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<dType, Type>* f, derivativeCol<dType>& D, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid,My_Int N_functions ,dType* params){ }

template<>
void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*b_cond) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_eom)(double *elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_b_cond)(double *elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid, grid<double>& Ygrid, My_Int N_functions, double* params);

template<>
void UpdateFunctions<double, Complex>( void (*eom) (Complex* elem, Cfunction<double, Complex> *F1 , grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*b_cond) (Complex* elem, Cfunction<double, Complex> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_eom)(Complex *elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
					 ,void (*l_b_cond)(Complex *elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<double, Complex>* f, derivativeCol<double>& D, grid<double>& Rgrid, grid<double>& Xgrid, grid<double>& Ygrid,My_Int N_functions,double* params);


template<>
void UpdateFunctions<float128, float128>( void (*eom) (float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*b_cond) (float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3,float128* params, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*l_eom)(float128 *elem, derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*l_b_cond)(float128 *elem,  derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<float128, float128>* f, derivativeCol<float128>& D, grid<float128>& Rgrid,grid<float128>& Xgrid, grid<float128>& Ygrid,My_Int N_functions, float128* params);

template<>
void UpdateFunctions<float128, ComplexQ>( void (*eom) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1 , grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*b_cond) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3,float128* params, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*l_eom)(ComplexQ *elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
										 ,void (*l_b_cond)(ComplexQ *elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, grid<float128>& grid2, grid<float128>& grid3, float128* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<float128, ComplexQ>* f, derivativeCol<float128>& D, grid<float128>& Rgrid, grid<float128>& Xgrid, grid<float128>& Ygrid,My_Int N_functions,float128* params);


template<>
void UpdateFunctions<mpreal, mpreal>( void (*eom) (mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k)
									 ,void (*b_cond) (mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3,mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k)
									 ,void (*l_eom)(mpreal *elem, derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
									 ,void (*l_b_cond)(mpreal *elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<mpreal, mpreal>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid,grid<mpreal>& Xgrid, grid<mpreal>& Ygrid,My_Int N_functions, mpreal* params);

template<>
void UpdateFunctions<mpreal, ComplexMP>( void (*eom) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1 , grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k)
										,void (*b_cond) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3,mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k)
										,void (*l_eom)(ComplexMP *elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
										,void (*l_b_cond)(ComplexMP *elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, grid<mpreal>& grid2, grid<mpreal>& grid3, mpreal* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<mpreal, ComplexMP>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid, grid<mpreal>& Xgrid, grid<mpreal>& Ygrid,My_Int N_functions,mpreal* params);

template<>
void UpdateFunctions<long double, long double>( void (*eom) (long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i, const My_Int& j, const My_Int& k)
											   ,void (*b_cond) (long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3,long double* params, const My_Int& i, const My_Int& j, const My_Int& k)
											   ,void (*l_eom)(long double *elem, derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
											   ,void (*l_b_cond)(long double *elem,  derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<long double, long double>* f, derivativeCol<long double>& D, grid<long double>& Rgrid,grid<long double>& Xgrid, grid<long double>& Ygrid,My_Int N_functions, long double* params);

template<>
void UpdateFunctions<long double, ComplexLD>( void (*eom) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1 , grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i, const My_Int& j, const My_Int& k)
											 ,void (*b_cond) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3,long double* params, const My_Int& i, const My_Int& j, const My_Int& k)
											 ,void (*l_eom)(ComplexLD *elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
											 ,void (*l_b_cond)(ComplexLD *elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, grid<long double>& grid2, grid<long double>& grid3, long double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<long double, ComplexLD>* f, derivativeCol<long double>& D, grid<long double>& Rgrid, grid<long double>& Xgrid, grid<long double>& Ygrid,My_Int N_functions,long double* params);

template<typename dType, typename Type>
dType MaxEOM(void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid, dType* params, const My_Int& i, const My_Int& j, const My_Int& k),My_Int N_eoms,Cfunction<dType, Type> *F1, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid, dType* params){
    
	dType w2=0;
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();
    My_Int start, end;
    My_Int TotalLength=(Np-2)*Nx*Ny;
    
    int numtasks;
    int rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
    My_Int div=TotalLength/numtasks;
    start=rank*div;
    
    if(rank==numtasks-1){
        end=TotalLength;
    }
    else{
        end=(rank+1)*div;
    }
    
	#pragma omp parallel
 	{Type *telem;
 	telem= new Type[N_eoms];
 	dType tw1=0, tw2=0;
    My_Int i;
    My_Int j;
    My_Int k;
	
 	#pragma omp for
 	for (My_Int ii=0; ii<end-start; ii++) {
        i=1+(ii+start)/(Nx*Ny);
        j=((ii+start) % (Nx*Ny))/Ny;
        k=((ii+start) % (Nx*Ny)) % Ny;

        eom(telem,F1, Rgrid, Xgrid, Ygrid, params, i, j,k);
        
        for (My_Int l=0;l<N_eoms; l++){
 				tw1=modulus<dType, Type>(telem[l]);
 				tw2=tw1>tw2?tw1:tw2;
            }}

	delete[] telem;
	
	#pragma omp critical
	{w2= tw2>w2?tw2:w2;}
	}
	
    dType *w_array;
    w_array=new dType[numtasks];
    
    MPI_Allgather(&w2, sizeof(dType),MPI_BYTE, w_array, sizeof(dType), MPI_BYTE,MPI_COMM_WORLD);
    
    for (int ii=0; ii<numtasks; ii++) {
        w2= w2>w_array[ii]?w2:w_array[ii];
    }
    
    delete[] w_array;
    
	return w2;
}

template<typename dType, typename Type>
dType MaxBCond(void (*b_cond) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid,dType* params, const My_Int& i, const My_Int& j, const My_Int& k), My_Int N_bconds,Cfunction<dType, Type> *F1,dType* params, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid, My_Int position){
	dType w2=0;
	My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();
    My_Int start, end;
    My_Int TotalLength=Nx*Ny;
    
    int numtasks;
    int rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    My_Int div=TotalLength/numtasks;
    start=rank*div;
    
    if(rank==numtasks-1){
        end=TotalLength;
    }
    else{
        end=(rank+1)*div;
    }
    
	
	#pragma omp parallel
	{Type *telem;
	telem= new Type[N_bconds];
	dType tw1=0, tw2=0;
    My_Int j;
    My_Int k;
	
	#pragma omp for
	for (My_Int ii=0; ii<end-start; ii++) {
        j=start/Ny;
        k=start % Ny;
		b_cond(telem, F1, Rgrid, Xgrid, Ygrid, params, position, j, k);
		for (My_Int l=0;l<N_bconds; l++){
			tw1= modulus<dType, Type>(telem[l]);
			tw2=tw1>tw2?tw1:tw2;
        }}
	delete [] telem;

    //note: "critical"
    //here it find the biggest w2 for each threat
    #pragma omp critical
    {w2= tw2>w2?tw2:w2;}
	
	}
    
    dType *w_array;
    w_array=new dType[numtasks];
    
    //recombines, and find the biggest out of all threats
    MPI_Allgather(&w2, sizeof(dType),MPI_BYTE, w_array, sizeof(dType), MPI_BYTE,MPI_COMM_WORLD);
    for (int ii=0; ii<numtasks; ii++) {
        w2= w2>w_array[ii]?w2:w_array[ii];
    }
    
    
    delete[] w_array;
    
	return w2;
}

template<typename dType, typename Type>
void CheckEoms(void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid, dType* params, const My_Int& i, const My_Int& j, const My_Int& k), My_Int N_eoms, Cfunction<dType, Type>*Fsol,Cfunction<dType, Type> *Ftest, grid<dType>& Rgrid, grid<dType>& Xgrid, grid<dType>& Ygrid, dType*params){

	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();

	#pragma omp parallel
	{Type *telem;
	telem=new Type[N_eoms];
	
        //Ftest is an array of Cfunctions
	#pragma omp for
	for (My_Int i=1;i<Np-1;i++){
		for (My_Int j=0;j<Nx;j++){
            for(My_Int k; j<Ny; k++){
			eom(telem,Fsol,Rgrid, Xgrid, Ygrid,params,i,j,k);
			for(My_Int l=0; l<N_eoms; l++){
				Ftest[l](i,j,k)=telem[l];//equations of motions for all fields at all points
            }
			}
		}
	}
	
	delete[] telem;}
}

#endif /* defined(__PDEs__NewtonMethod__) */
