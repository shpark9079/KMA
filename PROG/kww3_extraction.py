#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 08:50:20 2020

@author: shpark
"""

from scipy.io import netcdf
import numpy as np
import glob
import pandas as pd
import sys
import kma_module

args=sys.argv[1:]
base_dir=args[0]
infile=args[1]
srtdate=args[2]

base_dir='/data1/2020/KMA/SFR-01/1'

infile='kww3.2020021005.nc'

ncid=netcdf.netcdf_file(base_dir+'/'+infile,'r')

lon=ncid.variables['longitude'][:]
lat=ncid.variables['latitude'][:]
time=ncid.variables['time'][:]
hs=ncid.variables['WVHGT_surface']
wd=ncid.variables['WVDIR_surface'][:]
wp=ncid.variables['PWPER_surface'][:]
u=ncid.variables['UGRD_surface'][:]
v=ncid.variables['VGRD_surface'][:]
miss_v=hs._FillValue

col_name=['stnid','lat','lon','name']

yy='{0:04d}'.format(2020)
mon='{0:02}'.format(2)
dd='{0:02d}'.format(10)
hh='{0:02d}'.format(5)
mm='{0:02d}'.format(0)
srtdate=yy+mon+dd+hh+mm
for loc in glob.glob(base_dir+'/DAIN/*.dat'):
#for loc in glob.glob(base_dir+'/DAIN/sea_buoy_khoa.dat'):
    loc1=pd.read_csv(loc,header=None,sep='\t')
    loc1.columns=col_name
    for i in range(len(loc1['stnid'])):
#    for i in range(3,4):
        lon1=loc1['lon'][i]
        lat1=loc1['lat'][i]
        stn_id=str(loc1['stnid'][i])
        name=loc1['name'][i]
        loc_i,loc_j=kma_module.fine_loc(lon,lat,lon1,lat1,hs[0,:,:],miss_v)        
        hs1=hs[loc_i,loc_j]            
        outfile=base_dir+'/DASV/KWW3-'+name+'-'+stn_id+'-'+srtdate+'.dat'
        f=open(outfile,'w')
        f.write('YEAR MONTH DAY HOUR MINUTE HS\n')
        f.write('%s %s %s %s %s %3.4f\n' % (yy,mon,dd,hh,mm,hs1))
        f.close()
