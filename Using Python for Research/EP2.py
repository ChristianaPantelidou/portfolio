#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 16:46:37 2019

@author: christianapantelidou
"""

s=0;
f1=1;
f2=1;
while f2<=4000000:
    fnew=f1+f2;
    if fnew%2 == 0:
        s=s+fnew
    f1=f2
    f2=fnew
print(s)