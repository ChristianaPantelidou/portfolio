#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 21:55:22 2019

@author: christianapantelidou
"""
def IsPrime(n):
    rr=True
    for j in range(2,n):
        if (n% j) == 0:
            rr=False
            break
    return rr


c=0
n=2
while c<10001:
    if IsPrime(n)==True:
        c=c+1
        MaxPrime=n
        print(n)
    n=n+1