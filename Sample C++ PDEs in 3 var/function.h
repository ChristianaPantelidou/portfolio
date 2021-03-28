
//

#ifndef __PDEs__function__
#define __PDEs__function__

#include <iostream>
#include <complex>
#include "math.h"
#include "MyTypes.h"
#include "CustomFunctions.h"
#include "mpi.h"

// Grid class definitions
//My_Int is a long int
template<typename Type>
class grid{
private:
	My_Int TotPatches;
	My_Int TotLength;
	My_Int* PatchLeftBoundary;
	My_Int* PatchRightBoundary;
	Type* gridarray;
	My_Int* PatchOrder;
    My_Int* Uppers;
    My_Int* Lowers;
public:
	~grid();//destructor
	
	grid(My_Int* Ends, Type *TempGrid, My_Int* Orders, My_Int NPatches=1);   // This is the constructor
	
	My_Int OrderAtPoint(const My_Int& i);
	
	My_Int MaxOrder();
	
	My_Int LeftPatchBoundary(const My_Int& i);
	
	My_Int RightPatchBoundary(const My_Int& i);
    
    My_Int Lower(const My_Int& i);
    
    My_Int Upper(const My_Int& i);
    
    void NewtonLowerUpper(const My_Int& i, My_Int &Lower, My_Int &Upper);
	
	My_Int TotalLength();
	
	bool IsLeftBoundary(const My_Int& i);
	
	bool IsRightBoundary(const My_Int& i);
	
	Type operator[](const My_Int& i);
};

//Grid class implementation

//called this way: rgrid(Rdimensions, trgrid, ROrders,rlength);
//Rdimensions:dimension of each patch
//trgrid: is a dynamically allocated vector that contains the rgrid (read from file)
//ROrders: order of interpolation on each patch
//rlength: number of patches
template<typename Type>
grid<Type>::grid(My_Int* Ends, Type *grid, My_Int* Orders, My_Int NPatches){
    TotPatches=NPatches;
    //length=total length
    My_Int length(0);
    for (My_Int ii=0; ii<TotPatches; ++ii) {
        length+=Ends[ii];
    }
    TotLength=length;
    PatchLeftBoundary= new My_Int[length];
    PatchRightBoundary= new My_Int[length];
    PatchOrder= new My_Int[length];
    gridarray=new Type[length];
    Uppers=new My_Int[length];
    Lowers=new My_Int[length];
    
    //Ends=Rdimensions
    //eg End={2,3,1}, NPatches=3 Length=6
    //tempEnds={2,5,6}
    My_Int *tempEnds;
    tempEnds= new My_Int[TotPatches];
    tempEnds[0]=Ends[0];
    for (My_Int ii=1; ii<TotPatches; ++ii) {
        tempEnds[ii]=tempEnds[ii-1]+Ends[ii];
    }
    
    //For 0th patch only
    //eg, for ii=0,1 PatchLeftBoundary={0,0},PatchRightBoundary={1,1}
    for (My_Int ii=0; ii<Ends[0]; ++ii) {
        PatchLeftBoundary[ii]=0;
        PatchRightBoundary[ii]=tempEnds[0]-1;
        PatchOrder[ii]=Orders[0];
        gridarray[ii]=grid[ii];
    }
    
    //For the rest of the patches
    //eg, for 1st patch: ii=0,1,2, PatchLeftBoundary={2,2,2},PatchRightBoundary={4,4,4}
     //eg, for 2nd patch: ii=0, PatchLeftBoundary={5},PatchRightBoundary={5}
    for(My_Int jj=1; jj<TotPatches; ++jj){
        for (My_Int ii=tempEnds[jj-1]; ii<tempEnds[jj]; ++ii) {
            gridarray[ii]=grid[ii];
            PatchOrder[ii]=Orders[jj];
            PatchLeftBoundary[ii]=tempEnds[jj-1];//end point of previous patch+1
            PatchRightBoundary[ii]=tempEnds[jj]-1;//end point of current patch-1
        }
    }
    
    //length=total length of grid
    //PatchLeftBoundary={0,0,2,2,2,5}
    //PatchRightBoundary={1,1,4,4,4,5}
    for (My_Int ii=0; ii<length; ii++) {
        
        //eg, length=6, lowers={0-N0,1-N0,2-N1,3-N1,4-N1,5-N2}, uppers={0+N0,1+N0,2+N1,3+N1,4+N1,5+N2}
        Lowers[ii]= ii-PatchOrder[ii];
        Uppers[ii]=ii+PatchOrder[ii];
        
        if (Lowers[ii]< PatchLeftBoundary[ii]) {
            Lowers[ii]=PatchLeftBoundary[ii];
            Uppers[ii]= (PatchLeftBoundary[ii]+2*PatchOrder[ii]+1 > PatchRightBoundary[ii] ) ? PatchRightBoundary[ii] : PatchLeftBoundary[ii]+2*PatchOrder[ii]+1 ;
        }
        if (Uppers[ii]> PatchRightBoundary[ii] ) {
            Uppers[ii]=PatchRightBoundary[ii];
            Lowers[ii]= (PatchRightBoundary[ii]-2*PatchOrder[ii] -1 < PatchLeftBoundary[ii]) ? PatchLeftBoundary[ii] : PatchRightBoundary[ii]-2*PatchOrder[ii] -1 ;
        }
    }
    
    delete []tempEnds;
}

template<typename Type>
grid<Type>::~grid(){
    delete []PatchLeftBoundary;
    delete []PatchRightBoundary;
    delete []PatchOrder;
    delete []gridarray;
    delete []Uppers;
    delete []Lowers;
    PatchLeftBoundary=0;
    PatchRightBoundary=0;
    PatchOrder=0;
    Uppers=0;
    Lowers=0;
    gridarray=0;
}

template<typename Type>
My_Int grid<Type>::TotalLength(){
    return TotLength;
}

template<typename Type>
My_Int grid<Type>::LeftPatchBoundary(const My_Int& i){
    return PatchLeftBoundary[i];
}

template<typename Type>
My_Int grid<Type>::RightPatchBoundary(const My_Int& i){
    return PatchRightBoundary[i];
}

template<typename Type>
bool grid<Type>::IsLeftBoundary(const My_Int& i){
    return PatchLeftBoundary[i]==i;
}

template<typename Type>
bool grid<Type>::IsRightBoundary(const My_Int& i){
    return PatchRightBoundary[i]==i;
}

template<typename Type>
Type grid<Type>::operator[](const My_Int& i){
    return gridarray[i];
}

template<typename Type>
My_Int grid<Type>::OrderAtPoint(const My_Int& i){
    return PatchOrder[i];
}

template<typename Type>
My_Int grid<Type>::MaxOrder(){
    My_Int N(0);
    for(My_Int ii=0; ii<TotLength;++ii){
        N=(N>=PatchOrder[ii])?N:PatchOrder[ii];//if N>=PatchOrder[ii]], then N=PatchOrder[ii]]
    }
    return N;
}

template<typename Type>
My_Int grid<Type>::Lower(const My_Int& i){
    return Lowers[i];
}

template<typename Type>
My_Int grid<Type>::Upper(const My_Int& i){
    return Uppers[i];
}


template<typename Type>
void grid<Type>::NewtonLowerUpper(const My_Int& i1,My_Int& lower, My_Int& upper){
    
    lower= i1-(*this).OrderAtPoint(i1);
    upper= i1+(*this).OrderAtPoint(i1)+1;
    
    if (lower-(*this).OrderAtPoint(i1)-1 <= (*this).LeftPatchBoundary(i1)) {
        lower= (*this).LeftPatchBoundary(i1)-1;
        upper= (*this).Upper(i1)+1;
    }
    
    if (upper+(*this).OrderAtPoint(i1)+1 >= (*this).RightPatchBoundary(i1) ) {
        upper= (*this).RightPatchBoundary(i1)+2;
        lower= (*this).Lower(i1);
        if (lower==(*this).LeftPatchBoundary(i1) ){
            lower--;
        }
    }
    
    //if it is the first point on the patch, but on the first patch
    if (i1==(*this).LeftPatchBoundary(i1) && i1>0) {
        lower=(*this).Lower(i1-1);
        upper=(*this).Upper(i1);
    }
    
}


//Derivative class definitions

template<typename Type>
class derivative{
private:
	grid<Type>* Rgrid;
	grid<Type>* Xgrid;
    grid<Type>* Ygrid;
	Type *pr; // saves the differentiation matices for each derivative.
	Type *px;
	Type *py;
public:
	derivative(Type *Pr=0, Type *Px=0, Type *Py=0,grid<Type>* TRgrid=0, grid<Type>* TXgrid=0, grid<Type>* TYgrid=0);
	
	derivative(const derivative& D);
	
	
    Type operator()(const My_Int& i1,const  My_Int& i2, const  My_Int& i3,const My_Int& j1,const My_Int& j2, const My_Int& j3);
	
	derivative& operator=(const derivative& D);
	
    ~derivative();//destructor
    
	bool getpr(){return (pr==0);};//Am I taking a derivative in this direction or not?
	bool getpx(){return (px==0);};
    bool getpy(){return (py==0);};
	
	My_Int getNp(){return (*Rgrid).TotalLength();}
	My_Int getNx(){return (*Xgrid).TotalLength();}
    My_Int getNy(){return (*Ygrid).TotalLength();}
	
	My_Int Rorder(const My_Int& i){return (*Rgrid).OrderAtPoint(i);}
	
	My_Int Xorder(const My_Int& i){return (*Xgrid).OrderAtPoint(i);}
    
    My_Int Yorder(const My_Int& i){return (*Ygrid).OrderAtPoint(i);}
	
	My_Int RLeftBoundary(const My_Int& i){return (*Rgrid).LeftPatchBoundary(i);}
	
	My_Int XLeftBoundary(const My_Int& i){return (*Xgrid).LeftPatchBoundary(i);}
    
    My_Int YLeftBoundary(const My_Int& i){return (*Ygrid).LeftPatchBoundary(i);}
	
	My_Int RRightBoundary(const My_Int& i){return (*Rgrid).RightPatchBoundary(i);}
	
	My_Int XRightBoundary(const My_Int& i){return (*Xgrid).RightPatchBoundary(i);}
    
    My_Int YRightBoundary(const My_Int& i){return (*Ygrid).RightPatchBoundary(i);}
    
    My_Int XLower(const My_Int& i){return (*Xgrid).Lower(i);}
    
    My_Int XUpper(const My_Int& i){return (*Xgrid).Upper(i);}
    
    My_Int RLower(const My_Int& i){return (*Rgrid).Lower(i);}
    
    My_Int RUpper(const My_Int& i){return (*Rgrid).Upper(i);}
    
    My_Int YLower(const My_Int& i){return (*Ygrid).Lower(i);}
    
    My_Int YUpper(const My_Int& i){return (*Ygrid).Upper(i);}
};



//Derivative class implementation
//called this way:  derivative<MyType> d100(cheby_d1,idx,idy,&rgrid,&xgrid,&ygrid);
//if you do take a derivative with respect to a coordinate, then keep/fill in the appropriate structures
//USED FOR: to define the derivative class at the begining of the code
template<typename Type>
derivative<Type>::derivative(Type *Pr, Type *Px, Type *Py,grid<Type>* TRgrid, grid<Type>* TXgrid, grid<Type>* TYgrid){
    Rgrid=TRgrid;
    Xgrid=TXgrid;
    Ygrid=TYgrid;
    My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
    My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
    My_Int Ny=(Ygrid==0)?0:(*Ygrid).TotalLength();
    
    pr = new Type[Np*Np];
    px = new Type[Nx*Nx];
    py = new Type[Ny*Ny];
    
    
    if (Px==0){
        delete []px;
        px=0;
    }
    else{
        for (My_Int i=0; i<Nx*Nx; i++) {
            px[i]=Px[i];
        }
    }
    
    if (Py==0){
        delete []py;
        py=0;
    }
    else{
        for (My_Int i=0; i<Ny*Ny ; i++) {
            py[i]=Py[i];
        }
    }
    
    if (Pr==0){
        delete []pr;
        pr=0;
    }
    else{
        for (My_Int i=0; i<Np*Np ; i++) {
            pr[i]=Pr[i];
        }
    }
}

//overload constructor
template<typename Type>
derivative<Type>::derivative(const derivative<Type>& D){ //D is a derivative. ie, it comes with all its private and public variables
    Rgrid=D.Rgrid;
    Xgrid=D.Xgrid;
    Ygrid=D.Ygrid;
    My_Int Np=(*Rgrid).TotalLength();
    My_Int Nx=(*Xgrid).TotalLength();
    My_Int Ny=(*Ygrid).TotalLength();
    
    pr = new Type[Np*Np];
    px = new Type[Nx*Nx];
    py = new Type[Ny*Ny];
    
    
    if (D.px==0){
        delete []px;
        px=0;
    }
    else{
        for (My_Int i=0; i<Nx ; i++) {
            for (My_Int j=0; j<Nx; j++) {
                px[i*Nx+j]=D.px[i][j];
            }}
    }
    
    if (D.py==0){
        delete []py;
        py=0;
    }
    else{
        
        for (My_Int i=0; i<Ny ; i++) {
            for (My_Int j=0; j<Ny; j++) {
                py[i*Ny+j]=D.py[i][j];
            }}
    }
    
    if (D.pr==0){
        delete []pr;
        pr=0;
    }
    else{
        for (My_Int i=0; i<Np ; i++) {
            for (My_Int j=0; j<Np; j++) {
                pr[i*Np+j]=D.pr[i][j];
            }}
    }
    
}

template<typename Type>
Type derivative<Type>::operator()(const My_Int& i1,const My_Int& i2,const My_Int& i3,const My_Int& j1,const My_Int& j2,const My_Int& j3){
    My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
    My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
    My_Int Ny=(Ygrid==0)?0:(*Ygrid).TotalLength();
    Type tempvar= Type(1);// intialise tempvar=1
    
    if(pr!=0){
        tempvar*=pr[i1*Np+j1];//if I take derivative
    }
    else{
        tempvar*=(i1==j1?1:0);//if I don't take derivative: if on the diagonal, then 1, else 0
    }
    
    if(px!=0){
        tempvar*=px[i2*Nx+j2];
    }
    else{
        tempvar*=(i2==j2?1:0);
    }
    
    if(py!=0){
        tempvar*=py[i3*Ny+j3];
    }
    else{
        tempvar*=(i3==j3?1:0);
    }
    
    if(pr==0&&px==0&&py==0){
        std::cout << "oops" << std::endl;
    }
    
    return tempvar;//tempvar contains D_{i1j1}D_{i2j2}D_{i3j3}
    //to take the derivative of the function f at the point {1i1,i2,i3} you have to sum over the j's
}

//to set a derivative equal to another
//eg d100=d200
template<typename Type>
derivative<Type>& derivative<Type>::operator=(const derivative<Type>& D){
    if (this==&D)//if the addresses match, no need to change anything
        return *this;
    
    //if not, make this equal to D
    My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
    My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
    My_Int Ny=(Ygrid==0)?0:(*Ygrid).TotalLength();
    
    if (pr==0){}
    else{
        delete[] pr;
    }
    
    if (px==0){}
    else{
        delete[] px;
    }
    
    if (py==0){}
    else{
        delete[] py;
    }
    
    pr=0;
    px=0;
    py=0;
    
    Rgrid=D.Rgrid;
    Xgrid=D.Xgrid;
    Ygrid=D.Ygrid;
    
    Np=(*Rgrid).TotalLength();
    Nx=(*Xgrid).TotalLength();
    Ny=(*Ygrid).TotalLength();
    
    pr = new Type[Np*Np];
    px = new Type[Nx*Nx];
    py = new Type[Ny*Ny];
    
    if (D.px==0){
        delete []px;
        px=0;
    }
    else{
        for (My_Int i=0; i<Nx*Nx ; i++) {
            px[i]=D.px[i];
        }
    }
    
    if (D.py==0){
        delete []py;
        py=0;
    }
    else{
        for (My_Int i=0; i<Ny*Ny ; i++) {
            py[i]=D.py[i];
        }
    }
    
    if (D.pr==0){
        delete []pr;
        pr=0;
    }
    else{
        for (My_Int i=0; i<Np*Np ; i++) {
            pr[i]=D.pr[i];
        }
    }
    return *this;
}

template<typename Type>
derivative<Type>::~derivative(){
    My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
    My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
    My_Int Ny=(Ygrid==0)?0:(*Ygrid).TotalLength();
    
    if (pr==0){}
    else{
        delete[] pr;
    }
    
    if (px==0){}
    else{
        delete[] px;
    }
    
    if (py==0){}
    else{
        delete[] py;
    }
    
    pr=0;
    px=0;
    py=0;
    Rgrid=0;
    Xgrid=0;
    Ygrid=0;
}



//Derivative collection class definitions

template<typename Type>
class derivativeCol{
private:
	grid<Type>* Rgrid;
	grid<Type>* Xgrid;
    grid<Type>* Ygrid;
	My_Int N_ders;
	derivative<Type> *Ders;
public:
	derivativeCol(My_Int n,grid<Type>* TRgrid=0, grid<Type>* TXgrid=0, grid<Type>* TYgrid=0);
	
	~derivativeCol();
	
	derivative<Type>& operator[](const My_Int &i);
	
	My_Int order();
	
	My_Int Rorder(const My_Int& i){return (*Rgrid).OrderAtPoint(i);}
	
	My_Int Xorder(const My_Int& i){return (*Xgrid).OrderAtPoint(i);}
    
    My_Int Yorder(const My_Int& i){return (*Ygrid).OrderAtPoint(i);}
	
	My_Int RLeftBoundary(const My_Int& i){return (*Rgrid).LeftPatchBoundary(i);}
	
	My_Int XLeftBoundary(const My_Int& i){return (*Xgrid).LeftPatchBoundary(i);}
    
    My_Int YLeftBoundary(const My_Int& i){return (*Ygrid).LeftPatchBoundary(i);}
	
	My_Int RRightBoundary(const My_Int& i){return (*Rgrid).RightPatchBoundary(i);}
	
	My_Int XRightBoundary(const My_Int& i){return (*Xgrid).RightPatchBoundary(i);}
    
    My_Int YRightBoundary(const My_Int& i){return (*Ygrid).RightPatchBoundary(i);}
	
};

//Derivative collection class implementation

template<typename Type>
derivativeCol<Type>::derivativeCol(My_Int n,grid<Type>* TRgrid, grid<Type>* TXgrid, grid<Type>* TYgrid ){
    Rgrid=TRgrid;
    Xgrid=TXgrid;
    Ygrid=TYgrid;
    N_ders=n;
    Ders= new derivative<Type>[N_ders];
}

template<typename Type>
derivativeCol<Type>::~derivativeCol(){
    delete[] Ders;
    Ders=0;
    Rgrid=0;
    Xgrid=0;
    Ygrid=0;
}

template<typename Type>
derivative<Type>& derivativeCol<Type>::operator[](const My_Int &i){return Ders[i];}

template<typename Type>
My_Int derivativeCol<Type>::order(){
    return N_ders;
}


//Function class definitions

template <typename  dType, typename Type>
class function{
private:
	
    My_Int Np;
    My_Int Nx;
    My_Int Ny;
    Type *f;
	
public:
	//overloading the function
    
	function(My_Int rad_points=1, My_Int Xboundary_points=1, My_Int Yboundary_points=1);
	
    function(Type ***array,My_Int rad_points, My_Int boundaryX_points, My_Int boundaryY_points);
    
    function(Type *array,My_Int rad_points, My_Int boundaryX_points, My_Int boundaryY_points);

	function(Type ***array,grid<dType> *TRgrid, grid<dType> *TXgrid, grid<dType> *TYgrid);
	
	function(const function &f1);
    //////
    
    ~function();
    
    Type& operator()(const My_Int& i, const My_Int& j, const My_Int& k);
	
	function& operator=(const function &origf);
	
	function operator+(const function &f1);
	
	function operator-(const function &f1);
	
	template<typename dType1, typename Type1>
	friend function<dType1,Type1> operator*(const Type1& factor, function<dType1,Type1> & f);
	
	template<typename dType1, typename Type1>
	friend std::ostream& operator<<(std::ostream &out, function<dType1,Type1> &f1);
    
    Type diff(derivative<dType> &diff_operator,const My_Int& i, const My_Int& j, const My_Int& k);
	
	function diff(derivative<dType> &diff_operator);
	
};

// Function class implementation
template<typename dType, typename Type>
function<dType, Type>::function(My_Int rad_points, My_Int boundaryX_points, My_Int boundaryY_points){
    Np=rad_points;
    Nx=boundaryX_points;
    Ny=boundaryY_points;
    My_Int array_size=Np*Nx*Ny;
    f=new Type[array_size];

    for (My_Int i=0; i<array_size; i++) {
            f[i]=Type(0);//Initialise to zero, 1d array
            }
}

template<typename dType, typename Type>
function<dType, Type>::function(Type ***array,My_Int rad_points, My_Int boundaryX_points, My_Int boundaryY_points){
    Np=rad_points;
    Nx=boundaryX_points;
    Ny=boundaryY_points;
    My_Int array_size=Np*Nx*Ny;
    f=new Type[array_size];
    
    My_Int temp_int;
    
    for (My_Int i=0; i<Np; i++) {
        temp_int=i*Nx*Ny-Ny;
        for (My_Int j=0; j<Nx; j++) {
            temp_int+=Ny;
            for (My_Int k=0; k<Ny; k++){
                f[temp_int+k]=array[i][j][k];//Initialise to a 3d array
            }
        }
    }
}

template<typename dType, typename Type>
function<dType, Type>::function(Type *array,My_Int rad_points, My_Int boundaryX_points, My_Int boundaryY_points){
    Np=rad_points;
    Nx=boundaryX_points;
    Ny=boundaryY_points;
    My_Int array_size=Np*Nx*Ny;
    f=new Type[array_size];
    
    for (My_Int i=0; i<array_size; i++) {
                f[i]=array[i];//Initialise to a 1d array
            }

}

//called this way: function<MyType, MyType> testf0(test0,&rgrid,&xgrid,&ygrid);
template<typename dType, typename Type>
function<dType,Type>::function(Type ***array, grid<dType>* TRgrid, grid<dType>* TXgrid, grid<dType>* TYgrid){
    Np=(*TRgrid).TotalLength();
    Nx=(*TXgrid).TotalLength();
    Ny=(*TYgrid).TotalLength();
    My_Int array_size=Np*Nx*Ny;
    f=new Type[array_size];
    
    My_Int temp_int;
    
    for (My_Int i=0; i<Np; i++) {
        temp_int=i*Nx*Ny-Ny;
        for (My_Int j=0; j<Nx; j++) {
            temp_int+=Ny;
            for (My_Int k=0; k<Ny; k++){
                f[temp_int+k]=array[i][j][k];//Initialise to a 3d array
            }
        }
    }
}


//overload
template<typename dType, typename Type>
function<dType,Type>::function(const function<dType,Type> &f1){
	Np=f1.Np;
	Nx= f1.Nx;
    Ny= f1.Ny;
    My_Int array_size=Np*Nx*Ny;
    f=new Type[array_size];
    
    for (My_Int i=0; i<array_size; i++) {
                f[i]=f1.f[i];
            }
}

//destructor
template<typename dType, typename Type>
function<dType,Type>::~function(){
    delete []f;
    f=0;
}

template<typename dType, typename Type>
Type& function<dType,Type>::operator()(const My_Int& i, const My_Int& j, const My_Int& k){return f[(i*Nx+j)*Ny+k];};


template<typename dType, typename Type>
function<dType, Type> function<dType,Type>::operator+(const function<dType,Type> &f1){
	My_Int array_size=Np*Nx*Ny;
	Type* temp;
    temp=new Type[array_size];
	
    for (My_Int i=0; i<array_size; i++) {
            temp[i]=f[i]+f1.f[i];
	}
    
	function<dType, Type> tempf(temp,Np,Nx,Ny);//creates a new function and stores the sum in there
	
	delete [] temp;
	
	return tempf;
};

template<typename dType, typename Type>
function<dType,Type> function<dType,Type>::operator-(const function<dType,Type> &f1){
    My_Int array_size=Np*Nx*Ny;
    Type* temp;
    temp=new Type[array_size];
    
    for (My_Int i=0; i<array_size; i++) {
        temp[i]=f[i]-f1.f[i];
    }
    
    function<dType, Type> tempf(temp,Np,Nx,Ny);
    
    delete [] temp;
    
    return tempf;
};

template<typename dType1, typename Type1>
function<dType1, Type1> operator*(const Type1& x,function<dType1,Type1>& f){
	function<dType1, Type1> tempf(f.Np,f.Nx,f.Ny);
    My_Int array_size=tempf.Np*tempf.Nx*tempf.Ny;
    
	for (My_Int i=0; i<tempf.Np; ++i) {
			tempf.f[i]=x*f.f[i]; //creates a new function and stores the result
            }
    
	return tempf;
}

template<typename dType, typename Type>
function<dType, Type>& function<dType, Type>::operator=(const function<dType, Type> &origf){
	if (this==&origf)
		return *this;
	
    delete []f;
    f=0;
	
	Np=origf.Np;
	Nx=origf.Nx;
	Ny=origf.Ny;
    My_Int array_size=Np*Nx*Ny;
	f=new Type[array_size];
	
	for (My_Int i=0; i<array_size; i++){
			f[i]=(origf.f)[i];// why not f[i]=origf.f[i];???
            }
	return *this;
}

//prints f1 in a mathematica style
template<typename dType, typename Type>
std::ostream& operator<<(std::ostream &out,function<dType, Type> &f1){
	out << "{ ";
	for (My_Int i=0; i<f1.Np; i++) {
		out<< "{ ";
		for (My_Int j=0; j<f1.Nx; j++) {
                out<< "{ ";
            for (My_Int k=0; k<f1.Ny; k++) {
			out<< f1(i,j,k);
			if (k<f1.Ny-1) {out << ", ";}
            }
			if (j<f1.Nx-1) {out << "}, ";}
		}
		if(i<f1.Np-1){out<< "}}, ";}
		out << std::endl;}
	out<< "}}}" << std::endl;
	return out;
}

//returns the value of the derivative of function, eg d100, at point i,j,k
template<typename dType, typename Type>
Type function<dType, Type>::diff(derivative<dType> &diff_operator, const My_Int& i, const My_Int& j, const My_Int& k){
	
	Type Sum=Type(0);
    My_Int Rlower;
    My_Int Rupper;
    My_Int Xlower;
    My_Int Xupper;
    My_Int Ylower;
    My_Int Yupper;
    
    //bool getpr(){return (pr==0);}//not taking derivative
	if(diff_operator.getpr()){
        Rlower=i;
        Rupper=i;
	}
    //  My_Int XUpper(const My_Int& i){return (*Xgrid).Upper(i);}
    else{
        Rlower= diff_operator.RLower(i);
        Rupper= diff_operator.RUpper(i);
    }
    
    if(diff_operator.getpx()){
        Xlower=j;
        Xupper=j;
    }
    else{
        Xlower= diff_operator.XLower(j);
        Xupper= diff_operator.XUpper(j);
    }
    
    if(diff_operator.getpy()){
        Ylower=k;
        Yupper=k;
    }
    else{
        Ylower= diff_operator.YLower(k);
        Yupper= diff_operator.YUpper(k);
    }
    

    for(My_Int i1=Rlower; i1<=Rupper; i1++){
        for(My_Int i2=Xlower; i2<=Xupper; i2++){
            for(My_Int i3=Ylower; i3<=Yupper; i3++){
                Sum+=diff_operator(i,j,k,i1,i2,i3)* (*this)(i1,i2,i3);
            }
        }
    }
	
	return Sum;
};

//returns the derivative of the function on the whole grid
template<typename dType, typename Type>
function<dType, Type> function<dType, Type>::diff(derivative<dType> &diff_operator){
    My_Int TotalLength=Np*Nx*Ny;
    Type *temp_global;
    Type *temp_global2;
    My_Int start, end;
    
    int numtasks;
    int rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);//how many cores
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    My_Int div=TotalLength/numtasks;//how many points each core has//last one has a few more

    //each process has a different rank, start, end
    start=rank*div;
    
    if(rank==numtasks-1){
        end=TotalLength;
    }
    else{
        end=(rank+1)*div;
    }
    
    temp_global=new Type[end-start];
    temp_global2=new Type[TotalLength];
    
    int *lengths;
    lengths=new int[numtasks];
    
    int *starts;
    starts=new int[numtasks];
    
    ////????????////for MPI gather!!!
    for (int ii=0; ii<numtasks-1; ii++) {
        lengths[ii]=sizeof(Type)*div;//how many bits each of the cores has to store
        starts[ii]=sizeof(Type)*ii*div;//how many bits preceed
    }
    
    //last core treat separately cause has a few more points
    lengths[numtasks-1]= sizeof(Type)*(TotalLength- (numtasks-1)*div);
    starts[numtasks-1]= sizeof(Type)*(numtasks-1)*div;
    ////////????????//////
    
#pragma omp parallel
{
    My_Int i;
    My_Int j;
    My_Int k;
#pragma omp for
    for (My_Int ii=0; ii<end-start; ii++) {
        i=(ii+start)/(Nx*Ny);
        j=((ii+start) % (Nx*Ny))/Ny;
        k=((ii+start) % (Nx*Ny)) % Ny;
                temp_global[ii]=(*this).diff(diff_operator, i, j, k);
            }

}
    
//        std::cout << "Hi from " << rank << ", I did from " << start << " to " << end << std::endl;
    
    MPI_Allgatherv (temp_global,(end-start) * sizeof(Type) ,MPI_BYTE, temp_global2, lengths,starts, MPI_BYTE,MPI_COMM_WORLD);
    
    //fill in this function tempf with the array temp_global2
    function tempf(temp_global2, Np, Nx,Ny);
    
    delete[] starts;
    delete[] lengths;
    delete[] temp_global;
    delete[] temp_global2;
    
    return tempf;
};

//Cached function class definitions


template <typename  dType, typename Type>
class Cfunction{
private:
    My_Int nf;
    function<dType, Type> *fs;
public:
    
    Cfunction(My_Int Nf=0);
    
    Cfunction(function<dType,Type> &F,derivativeCol<dType> &D);
    
    ~Cfunction();
    
    template<typename dType1, typename Type1>
    friend std::ostream& operator<<(std::ostream &out, Cfunction<dType1,Type1> &f1);
    
    Type& operator()(const My_Int& i, const My_Int& j, const My_Int& k);
    
    Cfunction& operator=(const Cfunction& f1);
    
    Cfunction operator-(const Cfunction& f1);
    
    Cfunction operator+(const Cfunction& f1);
    
    template<typename dType1, typename Type1>
    friend Cfunction<dType1, Type1> operator*(const Type1& factor, Cfunction<dType1, Type1> &f);
    
    template<typename dType1, typename Type1>
    friend Cfunction<dType1, Type1> operator*(const dType1& factor, Cfunction<dType1, Type1> &f);
    
    Type diff(const My_Int& di, const My_Int& i, const My_Int& j, const My_Int& k);
    
    function<dType, Type>& diff(const My_Int& di);
    
    void update(derivativeCol<dType> &D);
    
    
};



//Cached function class implementation



////Overload
template<typename dType, typename Type>
Cfunction<dType, Type>::Cfunction(My_Int Nf){
	nf=Nf;
	if(nf==0){fs=0;}
	else{
		fs=new function<dType, Type>[Nf];
	}
}

//Given F, it completes the CFunction
template<typename dType, typename Type>
Cfunction<dType, Type>::Cfunction(function<dType, Type> &F,derivativeCol<dType> &D){
	nf=1+D.order();//to inclide the function itself as well
	fs=new function<dType, Type>[nf];
	fs[0]=F;
//	#pragma omp parallel for
	for (My_Int i=0; i<nf-1; ++i) {
		fs[i+1]=F.diff(D[i]);
	}
}

//given that the fs[0] is correct, it computes the rest
template<typename dType, typename Type>
void Cfunction<dType, Type>::update(derivativeCol<dType> &D){
    //	#pragma omp parallel for
    for (My_Int i=0; i<nf-1; ++i) {
        fs[i+1]=fs[0].diff(D[i]);
    }
}
//////



//destructor
template<typename dType, typename Type>
Cfunction<dType, Type>::~Cfunction(){
    if(fs!=0){delete [] fs;}
    fs=0;
}

//returns function at point i,j,k
template<typename dType, typename Type>
Type& Cfunction<dType, Type>::operator()(const My_Int& i, const My_Int& j, const My_Int& k){return fs[0](i,j,k);}

//returns the derivative di function at point i,j,k
template<typename dType, typename Type>
Type Cfunction<dType, Type>::diff(const My_Int& di, const My_Int& i, const My_Int& j, const My_Int& k){return fs[di](i,j,k);}

//returns the derivative di as a function, ie at all points
template<typename dType, typename Type>
function<dType, Type>& Cfunction<dType, Type>::diff(const My_Int& di){
	return fs[di];
}

template<typename dType, typename Type>
Cfunction<dType, Type>& Cfunction<dType, Type>::operator=(const Cfunction<dType,Type> &f1){
	if (fs!=0){delete [] fs;};
	nf=f1.nf;
	fs=new function<dType, Type>[nf];
	for (My_Int i=0; i<nf; ++i) {
		fs[i]=(f1.fs)[i];
	}
	return *this;
}

template<typename dType, typename Type>
Cfunction<dType, Type> Cfunction<dType, Type>::operator-(const Cfunction<dType, Type>& f1){
	Cfunction<dType, Type> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=fs[i]-f1.fs[i];
	}
	tempCf.nf= f1.nf ;
	return tempCf;
}

template<typename dType1, typename Type1>
Cfunction<dType1, Type1> operator*(const Type1& x, Cfunction<dType1, Type1>& f1){
	Cfunction<dType1, Type1> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=x*f1.fs[i];
	}
	
	return tempCf;
}

//multiply by scalar
template<typename dType1, typename Type1>
Cfunction<dType1, Type1> operator*(const dType1& x, Cfunction<dType1, Type1>& f1){
	Cfunction<dType1, Type1> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=Type1(x)*f1.fs[i];
	}
	
	return tempCf;
}

template<typename dType, typename Type>
std::ostream& operator<<(std::ostream &out, Cfunction<dType, Type> &f1){
	out << (f1.fs)[0];
	return out;
}

#endif /* defined(__PDEs__function__) */
