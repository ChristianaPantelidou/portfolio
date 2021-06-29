//
//  main.cpp
//  
//
//  Created by Christiana Pantelidou on 29/06/2021.
//

#include <stdio.h>
#include <iostream>
using namespace std;
#include "Vector3DStack.hpp"

int main()
{
    
    Vector3DStack v1(1,1,1);
    Vector3DStack v2(2,2,2);
    
//    cout<<v1*2<< std::endl;
//    cout<<v1*v2<< std::endl;
//    cout<<v1%v2<< std::endl;
//    cout<<v1.extractx()<< std::endl;
//    cout<<v2.norm()<< std::endl;
//    cout<<v2+v1<< std::endl;
    v2=v1.unitVector();
    cout<<v2<< std::endl;
    return 0;
}
