#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 21:50:15 2019

@author: christianapantelidou
"""

s=0
s2=0
for i in range(101):
    print(i)
    s=s+i
    s2=s2+i*i

print(s*s-s2)