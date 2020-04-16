#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:42:29 2020

@author: shpark
"""

import numpy as np
import pandas as pd
from astropy.time import Time

def ecmwf_ens_area(index,index1,index2,var,time_len,*param):
    
    ind=[601,602,603,604,605]
    ind_du=np.where(index==ind)
    
    if (len(ind_du[0]) !=0):
        ind1=index2
    else:
        ind1=index1
        
    dumy=np.where(ind1==index)
    min2=[]
    max2=[]
    v75_1=[]
    v25_1=[]
    median2=[]
    v50_1=[]
    for i in range(time_len):
          var3=[]          
          for t in range(len(var)):
              var1=var[t,:,:,:]
              var2=var1[:,dumy[0],dumy[1]]
              var2=var2[~var2.mask]
              var3.append(list(var2))
              
          min1=np.nanmin(var3)
          max1=np.nanmax(var3)
          v75=np.nanpercentile(var3,75)
          v25=np.nanpercentile(var3,25)
          v50=np.nanpercentile(var3,50)
          median1=np.nanmedian(var3)
          min2.append(min1)
          max2.append(max1)
          v75_1.append(v75)
          v50_1.append(v50)
          v25_1.append(v25)
          median2.append(median1)
          
    
    return   min2,max2,v75_1,v25_1,v50_1,median2   


def ecmwf_ens_seaway(loc_i,loc_j,var,time_len,*param):
    
    min2=[]
    max2=[]
    v75_1=[]
    v25_1=[]
    median2=[]
    v50_1=[]
    for i in range(time_len):
        var2=var[i,:,loc_i,loc_j]
            
        min1=np.nanmin(var2)
        max1=np.nanmax(var2)
        v75=np.nanpercentile(var2,75)
        v50=np.nanpercentile(var2,50)
        v25=np.nanpercentile(var2,25)
        median1=np.nanmedian(var2)
        min2.append(min1)
        max2.append(max1)
        v75_1.append(v75)
        v50_1.append(v50)
        v25_1.append(v25)
        median2.append(median1)
        
    return min2,max2,v75_1,v25_1,v50_1,median2     