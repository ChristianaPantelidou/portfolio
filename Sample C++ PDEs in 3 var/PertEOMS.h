
//

#ifndef __TimeDependence__TimeEvoEOMS__
#define __TimeDependence__TimeEvoEOMS__

#include <iostream>
#include "function.h"
#include "CustomFunctions.h"
#include "MyTypes.h"

template<typename dType, typename Type>
void equations(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void bconditions(Type* elem, Cfunction<dType, Type> *F1,dType* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<typename dType, typename Type>
void SecOrderEqs(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, const My_Int& i, const My_Int& j, const My_Int& k);







template<>
void equations<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1, grid<double>& grid1, double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1,double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void SecOrderEqs<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1, grid<double>& grid1, double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1,float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void SecOrderEqs<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, float128* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1,mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void SecOrderEqs<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, mpreal* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void equations<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void bconditions<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1,long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

template<>
void SecOrderEqs<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, long double* params, const My_Int& i, const My_Int& j, const My_Int& k);

#endif /* defined(__TimeDependence__TimeEvoEOMS__) */
