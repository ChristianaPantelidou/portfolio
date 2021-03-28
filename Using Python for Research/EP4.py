#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 19:19:29 2019

@author: christianapantelidou
"""

def IsPalindrom(n):
    result=False
    split=[i for i in str(n)]
    if split[0]==split[5] and split[1]==split[4] and split[2]==split[3]:
        result=True
    return result

set=[]
for i in range(999,100,-1):
    for j in range(999,100,-1):
        prod=i*j
        if prod>100000:
            if IsPalindrom(prod)==True:
                set.append(prod)
                break
max(set)