#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 19:38:29 2019

@author: christianapantelidou
"""

p=1*2*3*5*7*11*13*17*19;
found=False
while found==False:
    for i in range(1,20):
        if p%i!=0:
            break
        result=p
        found=True
    p=p+1
    
      