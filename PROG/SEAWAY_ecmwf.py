#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 17:34:14 2019

@author: shpark
"""

from scipy.io  import netcdf
from astropy.io import ascii
import sys
import numpy as np
from astropy.time import Time
import wind_dir_cal

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
seaway_index=dain+'/seaway_only_index.dat'

ecmwf_nc=netcdf.netcdf_file(daio+'/ecmwf_ense1.nc','r')
#ecmwf_nc=netcdf.netcdf_file(daio+'/'+netcdf_file,'r')
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


time=ecmwf_nc.variables['time']
st=time.units[12:]
st1=Time(st,format='iso',scale='utc')
time_len=len(time[:])
index=ascii.read(seaway_index)
for i in range(len(index['col1'])):
    loc_i=index['col4'][i]
    loc_j=index['col5'][i]
    hs3=hs2[:,loc_i,loc_j]
    wd3=wd2[:,loc_i,loc_j]
    wp3=wp2[:,loc_i,loc_j]
    wind_vel3=wind_vel2[:,loc_i,loc_j]
    wind_dir3=wind_dir2[:,loc_i,loc_j]
    out_file=dasv+'/SEAWAY_ecmwf_only_'+str(index['col1'][i])+'.dat'
    f=open(out_file,'w')
    f.write('Year Month Day Hour Hs WD WP WIND_vel WIND_dir\n')
    for j in range(len(time[:])):
		nt=Time(time[j]*3600,format='unix',scale='utc')
        nt1=st1.unix+nt.unix
        nt2=Time(nt1,format='unix',scale='utc')

        yy=nt2.iso[0:4]
        mm=nt2.iso[5:7]
        dd=nt2.iso[8:10]
        hh=nt2.iso[11:13]
        f.write('%s %s %s %s %3.4f %3.4f %3.4f %3.4f %3.4f\n' % (yy,mm,dd,hh,hs3[j],wd3[j],wp3[j],wind_vel3[j],wind_dir3[j]))
        
    f.close()
        
