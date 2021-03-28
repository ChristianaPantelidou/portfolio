
//

#include "NewtonMethod.h"

PetscErrorCode PCMGSetupViaCoarsen(PC pc,DM da_fine)
{
    PetscInt       nlevels,k,PETSC_UNUSED finest;
    DM             *da_list,*daclist;
    Mat            R;
    PetscErrorCode ierr;
    
    PetscFunctionBeginUser;
    nlevels = 1;
    PetscOptionsGetInt(NULL,NULL,"-levels",&nlevels,0);
    
    ierr = PetscMalloc1(nlevels,&da_list);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) da_list[k] = NULL;
    ierr = PetscMalloc1(nlevels,&daclist);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) daclist[k] = NULL;
    
    /* finest grid is nlevels - 1 */
    finest     = nlevels - 1;
    daclist[0] = da_fine;
    PetscObjectReference((PetscObject)da_fine);
    ierr = DMCoarsenHierarchy(da_fine,nlevels-1,&daclist[1]);CHKERRQ(ierr);
    for (k=0; k<nlevels; k++) {
        da_list[k] = daclist[nlevels-1-k];
//        DMDASetUniformCoordinates(da_list[k],0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
    }
    
    ierr = PCMGSetLevels(pc,nlevels,NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc,PC_MG_MULTIPLICATIVE);CHKERRQ(ierr);
    ierr = PCMGSetGalerkin(pc,PETSC_TRUE);CHKERRQ(ierr);
    
    for (k=1; k<nlevels; k++) {
        ierr = DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);CHKERRQ(ierr);
        ierr = PCMGSetInterpolation(pc,k,R);CHKERRQ(ierr);
        ierr = MatDestroy(&R);CHKERRQ(ierr);
    }
    
    /* tidy up */
    for (k=0; k<nlevels; k++) {
        ierr = DMDestroy(&da_list[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(da_list);CHKERRQ(ierr);
    ierr = PetscFree(daclist);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

  template<>
  void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*b_cond) (double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*l_eom)(double *elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k)
                                       ,void (*l_b_cond)(double *elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i1, const My_Int& j1, const My_Int& k1, const My_Int& i, const My_Int& j, const My_Int& k),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid, grid<double>& Ygrid,My_Int N_functions, double* params){
    

      
      int numtasks;
      int rank;
      
      double omp_start, omp_end, time, mtime;
      double omp_start_2, omp_end_2, time_2, mtime_2;
      
      MPI_Comm_size(PETSC_COMM_WORLD, &numtasks);
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      
    My_Int Np=Rgrid.TotalLength();
    My_Int Nx=Xgrid.TotalLength();
    My_Int Ny=Ygrid.TotalLength();
    My_Int nbdof = (Np-1)*Nx*Ny*N_functions; // number of degrees of freedom.
    My_Int N_lattice_sites=(Np-1)*Nx*Ny;// total, all processors
 	My_Int nvars = Ny*Nx+ (2*Rgrid.MaxOrder()+1)*(Nx+Ny)+N_functions*(2*Rgrid.MaxOrder()+1+Nx+Ny);
      
      std::vector<Eigen::Triplet<double, My_Int> > TriplVec;
      Eigen::Matrix<double, Eigen::Dynamic, 1> B(nbdof);//one for each line of the LinearOP
      
      TriplVec.reserve((long long)(nvars)* (long long)(nbdof)/(long long)numtasks);//each processor has its
      
      double tole=1.E-9;//Tolerance for the elements of the linear operator
      
      Vec            x, b, u, X;          /* approx solution, RHS, exact solution */
      VecScatter    Scatter;
      Mat            A;                /* linear system matrix */
      KSP            ksp;              /* linear solver context */
      PC             pc;               /* preconditioner context */
      PetscReal      norm,tol=1.e-11;  /* norm of solution error */
      PetscErrorCode ierr;
      PetscInt       n(1), ti, tj,rstart,rend,nlocal, its, sten_size;
      PetscInt       const *rstarts;
      PetscScalar    tvalue;
      IS ix,iX;
      
      DM    da;//3d matrix
      PetscInt  *i_map;//map from global to local coordinates
      DMDALocalInfo info;
      PetscInt  rs, xs, ys, re, xe, ye, rm, xm, ym;//s:starting, e:ending,m:distance
      AO    ao;//ordering context
      
      //Creates an object that will manage the communication of three-dimensional regular array data that is distributed across some processors. da: resulting matrix
      //Each entry has Nf values?
      DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX, Np-1,
                   Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,N_functions,1,NULL,NULL,NULL, &da);
      
      
      //Gets information about a given distributed array and this processors location in it
      //info: type DMDALocalInfo
      DMDAGetLocalInfo(da,&info);
      //Gets the application ordering context for a distributed array.
      //ao type AO: manages mapping
      DMDAGetAO(da,&ao);
      
      //(x,y,z)->(r,x,y)
      rs=info.xs; re=rs+ info.xm;
      xs=info.ys; xe=xs+ info.ym;
      ys=info.zs; ye=ys+ info.zm;
      
      //std::cout<< "rank"<<rank<<",re="<< re << ", rs="<< rs <<'\n';
      //std::cout<< "rank"<<rank<<",xe="<< xe << ", xs="<< xs <<'\n';
     // std::cout<< "rank"<<rank<<",ye="<< ye << ", rs="<< ys <<'\n';
      
      
      //Returns the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.
      //Output:
      //x,y,z	- the corner indices (where y and z are optional; these are used for 2D and 3D problems)
      //m,n,p	- widths in the corresponding directions (where n and p are optional; these are used for 2D and 3D problems)
      //this gives the same information as "info" above
      DMDAGetCorners(da, &rs,&xs,&ys,&rm,&xm,&ym);
      re=rs+rm; xe=xs+xm; ye=ys+ym;  //bottom right corenrs

      
      //number of lattice sites per processor
      N_lattice_sites= info.xm*info.ym*info.zm;
      
      
      rstarts= new PetscInt[numtasks+1];
      
      //Gets the integer value for a particular option in the database.????????
      ierr=PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
      
     //Creates other global vectors from da
      DMCreateGlobalVector(da, &x);
      DMCreateGlobalVector(da, &b);
      DMCreateGlobalVector(da, &u);
      
      i_map = new PetscInt[nbdof];//equal length as x,b,u//nbdof = (Np-1)*Nx*Ny*N_functions;--> total dof, not by each processor
      
      for (PetscInt ii=0; ii< nbdof; ++ii) {// ii global,flattened index
          i_map[ii]=ii;
      }
      
      //Maps a set of integers in the application-defined ordering to the PETSc ordering.
      //	ao	- the application ordering context
      //n	- the number of integers
      //ia	- the integers; these are replaced with their mapped value
      AOApplicationToPetsc(ao, nbdof, i_map);//da: 3d matrix with(Np-1)*Nx*Ny entries, each with length N_functions N//nbdof = (Np-1)*Nx*Ny*N_functions;

      
      
      VecGetOwnershipRange(x,&rstart,&rend);//rstart and rend refer to the flat vector
      VecGetOwnershipRanges(x,&rstarts);//rstarts contains the values of the rstart for each processor
      
//      std::cout << "Rank\t" << rank << "\t" << rstart << "\t" << rend-1 << std::endl;
      
      nlocal=rend-rstart;// length of the flat vector treated by each processor
      
      //Create this on each processor
      std::vector<My_Int> lindices;
      lindices.reserve(N_lattice_sites);//N_lattice_sites: number of lattice sites per processor
      
      for (My_Int kk=ys; kk<ye; ++kk) {
          for (My_Int jj=xs; jj<xe; ++jj){
              for (My_Int ii=rs; ii<re; ++ii) {
                  lindices.push_back(ii+(Np-1)*jj+(Np-1)*Nx*kk);//this is the set of flattened indices treated by each processor after being flattened
              }}}
      
      
//      std::cout << "Rank\t" <<rank  << "\tapp start\t" << i_map[N_functions*lindices[0]] << "\tapp end\t" << i_map[N_functions*lindices[lindices.size()-1]+N_functions-1] << std::endl;
      
      My_Int temp_int;
      
      
      omp_start = omp_get_wtime();
      
      //in the .h file. Constructs the linear operator. To do so uses precision "tole".-->TriplVec
      //Also evaluates the eom-->B
      ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, Ygrid, N_functions,params, lindices, i_map , tole );
      //TriplVec contains the column and row (global) of the Linear operator and the value, but only for the values in the particular processor

      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      //MPI down to the mother process only
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
            if (rank==0) {
      std::cout<<"Got linear operator in: "<< mtime <<'\n';
            }
      
      omp_start = omp_get_wtime();
      
      long long nz=TriplVec.size();//number of non-zeros in each processor
    
//      std::cout << "Rank: " << rank << "\tReserved:" << (long long)(nvars)* (long long)(nbdof)/(long long)numtasks << "\tNeeded: " << nz << std::endl;

      //here trying to count the non-zero elements: how?!
      My_Int *dnz;
      My_Int *onz;
      My_Int *totnz;
      My_Int *tdnz;
      My_Int *tonz;
      
      dnz=new My_Int[nbdof];//each is a 1d array! 1 for each line of the LinearOp
      onz=new My_Int[nbdof];
      totnz=new My_Int[nbdof];
      tdnz=new My_Int[nbdof];
      tonz=new My_Int[nbdof];
      
#pragma omp parallel for
for (My_Int ii=0; ii<nbdof; ++ii) {//initialises everything to zero
          dnz[ii]=0; onz[ii]=0;
          tdnz[ii]=0; tonz[ii]=0;
          totnz[ii]=0;
      }

      
My_Int *ldnz, *ltotnz;
      
#pragma omp parallel
{
    const int nthreads = omp_get_num_threads();//each processor is opening up threads
    const int ithread = omp_get_thread_num();//id number of each thread
    
#pragma omp single
{
    ldnz= new My_Int[nthreads*nbdof];//nbdof = (Np-1)*Nx*Ny*N_functions
    ltotnz= new My_Int[nthreads*nbdof];
}
    

#pragma omp for
    for (My_Int ii=0; ii<nthreads*nbdof; ++ii) {//initialises things
        ldnz[ii]=0;
        ltotnz[ii]=0;
    }
    
#pragma omp for
for (long long ii=0; ii<nz; ++ii) {
    if((rstart-1< TriplVec[ii].col() && TriplVec[ii].col()  < rend )){
        ldnz[TriplVec[ii].row()+ithread*nbdof]+=1;
    }
         ltotnz[TriplVec[ii].row()+ithread*nbdof]+=1;//can maximum fill in nlocal values
    
      }
    
#pragma omp for
for (My_Int ii=rstart; ii<rend; ++ii) {
            for (My_Int nn=0; nn<nthreads; nn++) {
            dnz[ii]+=ldnz[ii+nn*nbdof];
            totnz[ii]+=ltotnz[ii+nn*nbdof];
            }
        }
#pragma omp single
{
    delete[] ldnz;
    delete[] ltotnz;
}
    
}
      
      My_Int N_non_empty_rows(0);
      
#pragma omp parallel for reduction(+:N_non_empty_rows)
      for (My_Int ii=rstart; ii<rend; ii++) {
          if (totnz[ii]) {
          onz[ii]=totnz[ii]-dnz[ii];
            N_non_empty_rows++;
          }
      }
      
      My_Int **Col_Indices;
      PetscScalar **Vals;
      
      Col_Indices= new My_Int*[nbdof];
      Vals= new PetscScalar*[nbdof];
      
      My_Int temp_int1(0);
      
#pragma omp parallel private(temp_int1)
{
#pragma omp for
for (My_Int ii=0; ii<nbdof; ii++) {//goes through all the rows
          temp_int1=totnz[ii];//numbers of non-zeros on each row
          Col_Indices[ii]=0;
          Vals[ii]=0;
          if (temp_int1) {//if the line is not empty, go and open temp_int1slot in memory for the columns and for the values
              Col_Indices[ii]=new My_Int[temp_int1];
              Vals[ii]=new PetscScalar[temp_int1];
          }
      }
      }
      
      My_Int *temp_ints;
      
      temp_ints=new My_Int[nbdof];

#pragma omp parallel for
for (My_Int ii=0; ii<nbdof; ii++) {
          temp_ints[ii]=0;
      }
      
      omp_lock_t *lock;

      lock=new omp_lock_t[nlocal];
      
      
#pragma omp parallel for shared(lock)
      for (My_Int ii=0; ii<nlocal; ii++) {
          omp_init_lock(lock+ii);
      }
      

#pragma omp parallel for private(ti,tj, temp_int1) shared(lock, TriplVec, Col_Indices, Vals, temp_ints)
for(My_Int ii=0; ii<nz; ii++) {
          ti=TriplVec[ii].row();
          tj=TriplVec[ii].col();
        omp_set_lock(&(lock[ti-rstart]));//This subroutine forces the executing thread to wait until the specified lock is available. A thread is granted ownership of a lock when it becomes available.
        temp_int1=temp_ints[ ti ]++;//counter for each line, goes upto totnz[ti], less than nlocal
        omp_unset_lock(&(lock[ti-rstart]));
        Col_Indices[ ti ][ temp_int1 ] = tj;
        Vals[ ti ][ temp_int1 ] = TriplVec[ii].value(); //2d array
      }
      
      delete[] temp_ints;
      
      TriplVec.erase(TriplVec.begin(), TriplVec.end());
      
#pragma omp parallel for shared(lock)
      for (My_Int ii=0; ii<nlocal; ii++) {
          omp_destroy_lock(lock+ii);/////????????????????????????????????
      }
      
      delete []lock;
      
      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
      
      
      if (rank==0) {
          std::cout<<"Did loops in: "<< mtime  <<'\n';}
      
      MPI_Allreduce(dnz,tdnz,nbdof,MPI_LONG_LONG,MPI_SUM,PETSC_COMM_WORLD);//Combines values from all processes and distributes the result back to all processes
      MPI_Allreduce(onz,tonz,nbdof,MPI_LONG_LONG,MPI_SUM,PETSC_COMM_WORLD);
      
      My_Int av1(0),av2(0);
      
     MatCreateAIJ(PETSC_COMM_WORLD, nlocal,nlocal ,nbdof,nbdof, NULL, tdnz+rstart , NULL , tonz+rstart, &A);
      

      omp_start = omp_get_wtime();
      
      PetscInt i1,i2,i3;
      
      My_Int tii;
      
#pragma omp parallel private(tii, i3, i2, i1, ti, tvalue)
      {
#pragma omp for
          //at each point you get a 16*16 matrix for the linear Op and a 16 vector for b
          //insert a row in A and an element in B at the appropriate location
          //i1,i2,i3 are global indices
          //ti=local
          for(My_Int ii=0; ii<lindices.size(); ii++){
              tii=lindices[ii];//global indices treated by each process
              i3= tii/((Np-1)*Nx);
              i2= (tii%((Np-1)*Nx))/(Np-1);
              i1=1+(tii%((Np-1)*Nx))%(Np-1);
              for(My_Int l=0; l<N_functions; l++){
                  ti=i_map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l];//t1 is the row of the LInear operator matrix
              if (totnz[ti]) {
                  MatSetValues(A,1,&ti, totnz[ti] ,Col_Indices[ti],Vals[ti],INSERT_VALUES);
                  delete[] Col_Indices[ti];
                  delete[] Vals[ti];
              }
                  tvalue= B[ti];//ti's are the only non-zero on each processor
                  VecSetValues(b,1, &ti,&tvalue,INSERT_VALUES);
              }
          }
      }

      
      delete[] Col_Indices;
      delete[] Vals;
      
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);//reduce to a single process
      
      //put the global mtrices together
      MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
      VecAssemblyBegin(b);
      MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
      VecAssemblyEnd(b);
      
      omp_end = omp_get_wtime();
      
      time = omp_end - omp_start;
      
      MPI_Reduce(&time,&mtime,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);//reduce to to a single process
      
      if (rank==0) {
          std::cout<<"Assembled linear operator in: "<< mtime <<'\n';}

      //Creating the concept of the solver
      KSPCreate(PETSC_COMM_WORLD,&ksp);
      KSPGetPC(ksp,&pc);
      KSPSetOperators(ksp,A,A);//same matric for both KSP and PC
      KSPSetDM(ksp,da);
      KSPSetDMActive(ksp,PETSC_FALSE);//////Which method is used???????????????????
//      KSPSetType(ksp, KSPBCGS);
//      PCSetType(pc,PCBJACOBI);
      PCSetType(pc,PCMG);

      KSPSetTolerances(ksp,1.e-2,1.e-12,1.e+8,3000);//Sets the relative, absolute, divergence, and maximum iteration tolerances used by the default KSP convergence testers.
      
      PetscBool tf;
      
      PCSetFromOptions(pc);//pc=preconditioner
      KSPSetFromOptions(ksp);//setting the solver
      
      //Determines whether a PETSc object is of a particular type.
      PetscObjectTypeCompare((PetscObject)pc, PCMG,&tf);//tf=boolean
      
      if (tf) {
          PCMGSetupViaCoarsen(pc,da);} //not used!
      
      //Solve the actual system: Ax=b
      KSPSolve(ksp,b,x);
      
      KSPConvergedReason conv;
      KSPGetConvergedReason(ksp,&conv);

      KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
      
      //check result
      ierr = MatMult(A,x,u);//A x->u
      ierr = VecAXPY(u,-1.,b); //Computes u-b and stores it in u. Check how much is diviates from the original linear problem
      ierr = VecNorm(u,NORM_2,&norm);
      ierr = KSPGetIterationNumber(ksp,&its);//its=number of iterations
      
      PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);
      
      //Creates a standard, sequential array-style vector, X.
      VecCreateSeq(MPI_COMM_SELF,nbdof,&X);
      VecSetType(X,VECSEQ);//Equivalent with the above
      VecSetSizes(X,nbdof,nbdof);//Sets the local and global sizes, and checks to determine compatibility, global and local sizes respectively
      
      
      //Creates a data structure for an index set containing a list of evenly spaced integers.
      //nbdof	- the length of the locally owned portion of the index set
      //0	- the first element of the locally owned portion of the index set
      //1	- the change to the next index
      ISCreateStride(PETSC_COMM_WORLD ,nbdof,0,1, &iX);//iX is type IS
      
      //Creates a scatter context and executes it: x(result of linear solve)->X(empty), scattering fixed by the map ix->iX
      //x	- a vector that defines the shape (parallel data layout of the vector) of vectors from which we scatter
      //X	- a vector that defines the shape (parallel data layout of the vector) of vectors to which we scatter
      //ix	- the indices of x to scatter (if NULL scatters all values)
      //iX	- the indices of X to hold results (if NULL fills entire vector X)
      VecScatterCreate(x, NULL, X, iX, &Scatter);
      VecScatterBegin(Scatter,x,X,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(Scatter,x,X,INSERT_VALUES,SCATTER_FORWARD);
      
      
    //converts a petsc vector to an array
    PetscScalar *xx;
    VecGetArray(X,&xx);
      
      if(conv>0){
      for (My_Int l3=0; l3<N_functions; l3++) {
          for (My_Int i1=1; i1 < Np; i1++) {
              for (My_Int i2=0; i2<Nx; i2++) {
                  for (My_Int i3=0; i3<Ny; i3++){
                      f[l3](i1,i2,i3)=f[l3](i1,i2,i3)-xx[i_map[i3*N_functions*(Np-1)*Nx +i2*N_functions*(Np-1) +(i1-1)*N_functions +l3]];//updates the function
                  }}}}
      
      
      for (My_Int l3=0; l3<N_functions; ++l3) {
          f[l3].update(D);//Updates the derivatives
      }
      }
      else{
          std::cout << "Solver did not converge!" << std::endl;
          exit(1);
      }
      
      ierr = VecDestroy(&x); ierr = VecDestroy(&u);
      ierr = VecDestroy(&b); ierr = MatDestroy(&A);
      ierr = VecDestroy(&X); ISDestroy(&iX);
      ierr = KSPDestroy(&ksp);
//      AODestroy(&ao);
      DMDestroy(&da);
      VecScatterDestroy(&Scatter);
//      delete[] rstarts;
      delete[] dnz;
      delete[] onz;
      delete[] totnz;
      delete[] tdnz;
      delete[] tonz;
      delete[] i_map;

  }
