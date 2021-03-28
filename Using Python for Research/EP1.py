#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 16:39:05 2019

@author: christianapantelidou
"""

 s=0;
for i in range(1,1000): 
    if i%3 == 0 or i%5==0:
        s=s+i
print(s)