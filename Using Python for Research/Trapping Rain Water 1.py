#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 16:33:43 2019

@author: christianapantelidou
"""

      deriv=[]

        sign_set=[]
        for n in range(1,len(height)):
            if height[n-1]-height[n]<0:
                sign=-1
            elif height[n-1]-height[n]==0:
                sign=0
            else:
                sign=1
            deriv.append(height[n-1]-height[n])
            sign_set.append(sign)

        count=0
        min_set=[]
        left_set=[]
        right_set=[]
        h_right_set=[]
        h_left_set=[]
        n=0
        while n<len(deriv)-1:
            if deriv[n]<=0:
                n=n+1
                continue
            else:
                left=n
                minimum=n
                right=0
                for m in range(n+1, len(deriv)):
                    mm=m
                    if sign_set[m]==0 and m!=len(deriv)-1:
                        temp=[i for i, e in enumerate(sign_set[m+1:]) if e != 0]
                        if not temp:
                            mm=len(deriv)-1
                        else:
                            mm=(m+1)+temp[0]-1
                        temp=[]
                   # print(n,m,mm,sign_set)
                    if sign_set[m-1]-sign_set[mm]>0:
                        minimum=mm
                        right=len(deriv)
                        #print(n,m)
                        for s in range(m+1, len(deriv)):
                            t=s
                            #print(n,m,s)
                            if sign_set[s]==0 and s!=len(deriv)-1:
                               # print(sign_set[s+1:])
                                temp=[j for j, e in enumerate(sign_set[s+1:]
        ) if e != 0]
                               # print(temp[0])
                                if not temp:

                                    t=len(deriv)-1
                                else:
                                    t=(s+1)+temp[0]
                                #print(n,minimum,s,temp)
                                temp=[]
                            #print(n,minimum,s,t)
                            if sign_set[s-1]-sign_set[t]<0:
                                right=t
                                break
                        break
                if right==0:
                    break 
                left_set.append(left)
                right_set.append(right)
                h_left_set.append(height[left])
                h_right_set.append(height[right])
                min_set.append(minimum) 

                if right==0:
                    break
                else:
                    n=right
        #print(left_set)
       # print(right_set)
       # print(min_set)

        min_set=[]
        i=0
        while i<len(left_set):
            #print('here',i,h_left_set[i],max(h_right_set[i:]))
            if h_left_set[i]>max(h_right_set[i:]):
                index1=h_right_set[i:].index(max(h_right_set[i:])) #index with respect to i, this is the index in the lima
                   # print(i+index1)
                min_set.append((left_set[i],right_set[i+index1]))
                i=i+index1+1
            else: #find the first that is bigger or equal
                j=i;
                while j<len(left_set): #find the first right that is bigger
                   if h_right_set[j]>=h_left_set[i]:
                        min_set.append((left_set[i],right_set[j]))
                        i=j+1
                   # print('here2',i,j)
                        break
                   else:
                        j=j+1
               # print('exite loop')

        #print(min_set)            


        for (left,right) in min_set:
            threshold=min(height[left],height[right])
            for i in range(left,right):
               # print(left,right,i)
                if height[i]<threshold:
                    count=count+(threshold-height[i])