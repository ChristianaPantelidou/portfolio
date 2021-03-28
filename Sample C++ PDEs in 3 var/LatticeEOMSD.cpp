
//

#include "LatticeEOMS.h"


///The mapping between mathematica quantities and c quantities is done in the \\include EOMInFUncDefs line included in this code
//Here the functions are defined
#include "SourceFiles/EOMFuncDefs.c"

//returns the value of the eoms at a point i,j,k.. stored in elem
template<>
void equations<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,  double* params, const My_Int& i, const My_Int& j, const My_Int& k){
    
//call the functions "equationi" for i=1..16 one at a time to fill in the elem
#include "SourceFiles/EOMTotBody.c"
    
}

//correct elem for the boundary conditions
template<>
void bconditions<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3,double* params,  const My_Int& i, const My_Int& j, const My_Int& k){
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

#include "SourceFiles/BC.c"
}

template<>
void ksi2<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k){
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

#include "SourceFiles/Ksi2.c"
}

template<>
void psi<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, grid<double>& grid2, grid<double>& grid3, double* params, const My_Int& i, const My_Int& j, const My_Int& k){
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
    
#include "SourceFiles/Psi.c"
}
