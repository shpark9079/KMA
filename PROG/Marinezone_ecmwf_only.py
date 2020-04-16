#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:26:13 2019

@author: shpark
"""

from scipy.io  import netcdf
from astropy.io import ascii
import sys
import numpy as np
from astropy.time import Time
import numpy.ma as ma
import static_test
import re

#args=sys.argv[1:]
#daio=args[0]
#dain=args[1]
#dasv=args[2]
#srtdate=args[3]
#netcdf_file=args[4]

base_dir='/data1/2020/KMA/SFR-07'
daio=base_dir+'/DAIO'
dain=base_dir+'/DAIN'
dasv=base_dir+'/DASV'
srtdate='2020040800'
#area_index=base_dir+'/area_index.dat'
ecmwf_nc=netcdf.netcdf_file(daio+'/ecmwf_ense1.nc','r')
marine_zone=dain+'/marinezone.dat'

rww3_nc=netcdf.netcdf_file(netcdf_file,'r')
hs=ecmwf_nc.variables['swh']
wd=ecmwf_nc.variables['mwd']
wp=ecmwf_nc.variables['mwp']
wind_vel=ecmwf_nc.variables['wind']
wind_dir=ecmwf_nc.variables['dwi']
fill_v=hs._FillValue
hs1=ma.masked_where(hs[:,:,:,:]==hs._FillValue,hs[:,:,:,:])
wd1=ma.masked_where(wd[:,:,:,:]==wd._FillValue,wd[:,:,:,:])
wp1=ma.masked_where(wp[:,:,:,:]==wp._FillValue,wp[:,:,:,:])
wind_vel1=ma.masked_where(wind_vel[:,:,:,:]==wind_vel._FillValue,wind_vel[:,:,:,:])
wind_dir1=ma.masked_where(wind_dir[:,:,:,:]==wind_dir._FillValue,wind_dir[:,:,:,:])
hs2=(hs1*hs.scale_factor)+hs.add_offset
wd2=(wd1*wd.scale_factor)+wd.add_offset
wp2=(wp1*wp.scale_factor)+wp.add_offset
wind_vel2=(wind_vel1*wind_vel.scale_factor)+wind_vel.add_offset
wind_dir2=(wind_dir1*wind_dir.scale_factor)+wind_dir.add_offset
date=srtdate
marinezone=ascii.read(marine_zone)
index=marinezone
data=static_test.rww3_marinezone(dasv,marinezone,lon,lat,hs,wd,wp,u,v,time)
