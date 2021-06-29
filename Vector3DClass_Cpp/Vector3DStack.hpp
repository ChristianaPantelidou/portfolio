//
//  Vector3DStack.hpp
//  
//
//  Created by Christiana Pantelidou on 29/06/2021.
//

#ifndef Vector3DStack_hpp
#define Vector3DStack_hpp

#include <stdio.h>
#include <iostream>

#endif /* Vector3DStack_hpp */

class Vector3DStack{
    private:
        double x;
        double y;
        double z;
    
    public:
        Vector3DStack(const double x1,const double y1,const double z1);
        Vector3DStack(const Vector3DStack & V);
        ~Vector3DStack();
        double extractx() const;
        double extracty() const;
        double extractz() const;
        void setx(const double a);
        void sety(const double a);
        void setz(const double a);
        double norm() const;
        Vector3DStack operator=(const Vector3DStack  &V);
        Vector3DStack operator*(const double a) const;
        Vector3DStack operator/(const double a) const;
        Vector3DStack operator+(const Vector3DStack  &V) const;
        Vector3DStack operator-(const  Vector3DStack  &V) const;
        double operator*(const  Vector3DStack  &V) const;
        Vector3DStack operator%(const  Vector3DStack  &V) const;
        Vector3DStack unitVector() const;
        Vector3DStack getOrthogonalUnitVector(const Vector3DStack &V) const;
        friend std::ostream& operator<<(std::ostream &out,const Vector3DStack &V);
        bool operator==(const Vector3DStack &V) const;
        bool operator!=(const Vector3DStack &V) const;

};


