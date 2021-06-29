//
//  Vector3DStack.cpp
//  
//
//  Created by Christiana Pantelidou on 29/06/2021.
//

#include "Vector3DStack.hpp"
#include <stdio.h>
#include <iostream>
#include <math.h>

Vector3DStack::Vector3DStack(const double x1,const double y1,const double z1){//divide by scalar
    x=x1;
    y=y1;
    z=z1;
}

Vector3DStack::Vector3DStack(const Vector3DStack & V){//divide by scalar
    x=V.x;
    y=V.y;
    z=V.z;
}

Vector3DStack::~Vector3DStack(){//divide by scalar
    x=0;
    y=0;
    z=0;
}

double Vector3DStack::extractx() const{
    return x;
}
double Vector3DStack::extracty() const{
    return y;
}
double Vector3DStack::extractz() const{
    return z;
}
void Vector3DStack::setx(const double a){
    x=a;
}
void Vector3DStack::sety(const double a){
    y=a;
}
void Vector3DStack::setz(const double a){
    z=a;
}

Vector3DStack Vector3DStack::operator=(const Vector3DStack  &V){
    if (this==&V)
        return (*this);
    
    x=V.x;
    y=V.y;
    z=V.z;
    return (*this);
}

double Vector3DStack::norm() const{//find the norm
    return sqrt(x*x+y*y+z*z);
}
Vector3DStack Vector3DStack::operator*(const double a) const{// multiply by scalar
    Vector3DStack  Vnew(a*x, a*y, a*z);
    return Vnew;
}
Vector3DStack Vector3DStack::operator/(const double a) const{//divide by scalar
    if (a==0)
        throw std::runtime_error("Division by zero");
    Vector3DStack  Vnew(x/a, y/a, z/a);
    return Vnew;
}

Vector3DStack Vector3DStack::operator+(const  Vector3DStack & V) const{//add vectors
    Vector3DStack  Vnew(x+V.x, y+V.y, z+V.z);
    return Vnew;
}
Vector3DStack Vector3DStack::operator-(const  Vector3DStack& V) const{//subtract vectors
    Vector3DStack  Vnew(x-V.x, y-V.y, z-V.z);
    return Vnew;
}

double Vector3DStack::operator*(const Vector3DStack &  V) const{
    return (x)*(V.x)+(y)*(V.y)+(z)*(V.z);
}

Vector3DStack Vector3DStack::operator%(const Vector3DStack  & V) const{
    Vector3DStack Vnew(y*(V.z)-z*(V.y),z*(V.x)-x*(V.z),x*(V.y)-y*(V.x));
    return Vnew;
}

Vector3DStack Vector3DStack::unitVector() const{//subtract vectors
    double a=(*this).norm();
    Vector3DStack Vnew(x/a,y/a,z/a);
    return Vnew;
}

Vector3DStack Vector3DStack::getOrthogonalUnitVector(const Vector3DStack &V) const{
  Vector3DStack Vnew = (*this) % V;
  return Vnew.unitVector();
}

std::ostream& operator<<(std::ostream &out,const Vector3DStack &V){//print
    out<< "{"<<V.x<<","<<V.y<<","<<V.z<<"}"<< std::endl;
    return out;
}

bool Vector3DStack::operator==(const Vector3DStack &V) const
{
  return (x == V.x && y == V.y && z == V.z);
}

bool Vector3DStack::operator!=(const Vector3DStack &V) const
{
    return !operator==(V);
}

