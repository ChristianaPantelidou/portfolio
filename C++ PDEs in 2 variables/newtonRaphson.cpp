
#include "newtonRaphson.hpp"
#include <time.h>

void newtonRaphson(double * x,cFunction* f,double* xgrid,double*ygrid, derivativeClass & diffOperator,int Nf,int Nx,int Ny,double* c){
    double step;
    double *null = ( double * ) NULL;
    void *Numeric;
    int status;
    void *Symbolic;
    double* b= new double[Nf*Nx*(Ny-1)]();
    long * Ap=new long[Nf*Nx*(Ny-1)+1];
    long * Ai=new long[Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2)];
    double * Ax=new double[Nf*Nx*(Ny-1)*Nf*(Nx+Ny-2)];
    
    
    totalEqns(b,Nf,Nx,Ny,f,xgrid,ygrid,c);
    writeInFile("testingeom.bin",b,Nf*Nx*(Ny-1));
    
   
   variationEom(Ap,Ai,Ax,f,xgrid,ygrid,diffOperator,Nf,Nx,Ny,c);

  /*
   Eigen::VectorXd bb(Nf*Nx*(Ny-1)),xx(Nf*Nx*(Ny-1));
   for (int i=0;i<Nf*Nx*(Ny-1);i++){
   bb(i)=b[i];
   //std::cout<<b[i];
   }
   delete []b;
   
   A.makeCompressed();

   
  //Eigen::SparseLU<Eigen::SparseMatrix<double> >  solver;
    Eigen::UmfPackLU<Eigen::SparseMatrix<double> >  solver;
     solver.analyzePattern(A);
    // Compute the ordering permutation vector from the structural pattern of A
     solver.compute(A);
     solver.factorize(A);
     xx = solver.solve(bb);
    
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> >  solver;
   // solver.compute(A);
  //  xx = solver.solve(bb);
   
   
   for (int i=0;i<Nf*Nx*(Ny-1);i++){
   x[i]=xx(i);
   //std::cout<<x[i];
   }

   */
    
    
    //  Carry out the symbolic factorization.
    //
    status = umfpack_dl_symbolic (Nf*Nx*(Ny-1),Nf*Nx*(Ny-1), Ap, Ai, Ax, &Symbolic, null, null );
    //
    //  Use the symbolic factorization to carry out the numeric factorization.
    //
    status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null );
    //
    //  Free the memory associated with the symbolic factorization.
    //
    umfpack_dl_free_symbolic ( &Symbolic );
    //
    //  Solve the linear system.
    //
    status = umfpack_dl_solve ( UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null );
    //
    //  Free the memory associated with the numeric factorization.
    //
    umfpack_dl_free_numeric ( &Numeric );

    

  
      // t = clock() - t;
   // std::cout<<((double)t)/CLOCKS_PER_SEC;
    
    writeInFile("step.bin",x,Nf*Nx*(Ny-1));
}

        
