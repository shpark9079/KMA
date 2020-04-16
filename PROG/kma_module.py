#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:07:40 2020

@author: shpark
"""
import numpy as np

def fine_loc(model_lon,model_lat,lon,lat,var,*param):
    lon_dis=model_lon[1]-model_lon[0]
    lat_dis=model_lat[1]-model_lat[0]
    dumy_lon=np.where((model_lon > lon-lon_dis) & (model_lon < lon+lon_dis))
    dumy_lat=np.where((model_lat > lat-lat_dis) & (model_lat < lat+lat_dis))    
#    print(dumy_lon[0])
#    print(dumy_lat[0])
    lat2=[];lon2=[]
    for i in range(len(dumy_lon[0])):
      for j in range(len(dumy_lat[0])):
        if (var[dumy_lat[0][j],dumy_lon[0][i]] > 0):
            lon2.append(dumy_lon[0][i])
            lat2.append(dumy_lat[0][j])    
    lon1=model_lon[list(set(lon2))]
    lat1=model_lat[list(set(lat2))]
    
    dis_lon=lon1-lon
    dis_lat=lat1-lat
    dumy_lon=np.where(min(dis_lon)==dis_lon)
    dumy_lat=np.where(min(dis_lat)==dis_lat)

    lon_loc=np.where(model_lon==lon1[dumy_lon])
    lat_loc=np.where(model_lat==lat1[dumy_lat])
    if (var[lat_loc[0],lon_loc[0]] <100):
        return lat_loc[0],lon_loc[0]
    else:
        for i in range(0,2):
            for j in range(0,2):
                if (var[lat_loc[0]+i,lon_loc[0]+j] < 100):
                    return lat_loc[0]+i,lon_loc[0]+j
    