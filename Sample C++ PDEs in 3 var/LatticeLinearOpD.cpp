

#include "LatticeLinearOp.h"

//using namespce mpfr;


///definitions of 16*16 functions
//Includes the mapping
#include "SourceFiles/LOFuncDefs.c"



//fills in elem[1..256] for just boundary conditions in y and in between patches
void LOBulkLinearOp(double *elem, derivativeCol<double>& \
                            D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const \
                            My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
                            j1, const My_Int& k1){
#include "SourceFiles/LOInFuncDefs.c"//mapping
#include "SourceFiles/LO.c"
}

//gives elem[1..256] for the mixed derivatives
void LO2BulkLinearOp(double *elem, derivativeCol<double>& \
                       D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const \
                       My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
                       j1, const My_Int& k1){
#include "SourceFiles/LOInFuncDefs.c"//mapping
#include "SourceFiles/LO2.c"
}

//gives elem[1..256] for the variation of the rest of the terms in the equations of motion
//does not have the r=1 boundary conditions, ie BULK
//(i,j,k) is the point we are looking at!----> first point in the call!
template<>
void BulkLinearOp<double, double>(double *elem, derivativeCol<double>& D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){

    
    bool II(i==i1), JJ(j==j1), KK(k==k1);
    
    for(int ii=0; ii<256; ii++){
        elem[ii]=0;}
    
if((grid1.IsRightBoundary(i)) || (grid1.IsLeftBoundary(i)) || (k==0) || (k==grid3.TotalLength()-1)){
    
LOBulkLinearOp(elem, D, F1, grid1, grid2, grid3, params, i, j, k, i1, j1, k1);
}
else{
    
LO2BulkLinearOp(elem, D, F1, grid1, grid2, grid3, params, i, j, k, i1, j1, k1);
    
    if (II && JJ) {
//D001LOBulkLinearOp(elem, D, F1, grid1, params, i, j, k, i1, j1, k1);
        #include "SourceFiles/D001LOBody.c"
    }
    if (II && KK) {
//D010LOBulkLinearOp(elem, D, F1, grid1, params, i, j, k, i1, j1, k1);
        #include "SourceFiles/D010LOBody.c"
    }
    if (JJ && KK) {
//D100LOBulkLinearOp(elem, D, F1, grid1, params, i, j, k, i1, j1, k1);
        #include "SourceFiles/D100LOBody.c"
    }
    if (II && JJ && KK) {
//LocalLOBulkLinearOp(elem, D, F1, grid1, params, i, j, k, i1, j1, k1);
        #include "SourceFiles/LocalLOBody.c"
    }
    
}

}

//gives a elem[16*16] for the boundary
template<>
void BoundaryLinearOp<double, double>(double *elem, derivativeCol<double>& D,Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params,  const My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& j1, const My_Int& k1){
    Cfunction<double,double> &Qtt=F1[0];
    Cfunction<double,double> &Qrr=F1[1];
    Cfunction<double,double> &Q=F1[2];
    Cfunction<double,double> &B=F1[3];
    Cfunction<double,double> &R=F1[4];
    Cfunction<double,double> &Qr1=F1[5];
    Cfunction<double,double> &Qr2=F1[6];
    Cfunction<double,double> &Qt1=F1[7];
    Cfunction<double,double> &Qt2=F1[8];
    Cfunction<double,double> &Qtr=F1[9];
    Cfunction<double,double> &a0=F1[10];
    Cfunction<double,double> &a1=F1[11];
    Cfunction<double,double> &a2=F1[12];
    Cfunction<double,double> &ar=F1[13];
    Cfunction<double,double> &h1=F1[14];
    Cfunction<double,double> &h2=F1[15];
    grid<double>& rgrid=grid1;
    grid<double>& xgrid=grid2;
    grid<double>& ygrid=grid3;
    typedef double Type;
    double mu=params[0];
    double rh=params[1];
    double H=params[2];
    double R0=params[3];
    double L1=params[4];
    double L2=params[5];
    double S=params[6];
    double q=params[7];
    My_Int d100(1),d200(2),d110(3),d010(4),d020(5),d001(6),d002(7),d101(8),d011(9);
    
    for(int ii=0; ii<256; ii++){
        elem[ii]=0;}
	
#include "SourceFiles/BCLO.c"
	
}



//void D100LOBulkLinearOp(double *elem, derivativeCol<double>& \
//                        D,Cfunction<double, double> *F1, grid<double>& grid1,  double* params, const \
//                        My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
//                        j1, const My_Int& k1){
//#include "SourceFiles/LOInFuncDefs.c"
//#include "SourceFiles/D100LOFiles.c"
//}
//
//void D010LOBulkLinearOp(double *elem, derivativeCol<double>& \
//                           D,Cfunction<double, double> *F1, grid<double>& grid1,  double* params, const \
//                           My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
//                           j1, const My_Int& k1){
//#include "SourceFiles/LOInFuncDefs.c"
//#include "SourceFiles/D010LOFiles.c"
//}
//
//void D001LOBulkLinearOp(double *elem, derivativeCol<double>& \
//                           D,Cfunction<double, double> *F1, grid<double>& grid1,  double* params, const \
//                           My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
//                           j1, const My_Int& k1){
//#include "SourceFiles/LOInFuncDefs.c"
//#include "SourceFiles/D001LOFiles.c"
//}

//void LocalLOBulkLinearOp(double *elem, derivativeCol<double>& \
//                           D,Cfunction<double, double> *F1, grid<double>& grid1,  double* params, const \
//                           My_Int& i, const My_Int& j, const My_Int& k, const My_Int& i1, const My_Int& \
//                           j1, const My_Int& k1){
//#include "SourceFiles/LOInFuncDefs.c"
//#include "SourceFiles/LocalLOFiles.c"
//}
