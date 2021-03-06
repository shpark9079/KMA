#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:04:43 2020

@author: shpark
"""

from scipy.io  import netcdf
from astropy.io import ascii
import sys
import numpy as np
from astropy.time import Time
import numpy.ma as ma
import sub_ecmwf
#args=sys.argv[1:]
#daio=args[0]
#dain=args[1]
#dasv=args[2]
#srtdate=args[3]
#num=args[4]
#netcdf_file=args[5]


base_dir='/data1/2020/KMA/SFR-07'
daio=base_dir+'/DAIO'
dain=base_dir+'/DAIN'
dasv=base_dir+'/DASV'
srtdate='2020040800'
#area_index=base_dir+'/area_index.dat'
area_index=dain+'/area_index2.dat'

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
area=ascii.read(area_index)
ncid1=netcdf.netcdf_file(dain+'/ecmwf_ens_index1.nc')
ncid2=netcdf.netcdf_file(dain+'/ecmwf_ens_index2.nc')

index1=ncid1.variables['index'][:,:]
index2=ncid2.variables['index'][:,:]

for i in range(len(area['col1'])):
#for i in range(0,1):
    hs_mi,hs_ma,hs_75,hs_25,hs_50,hs_median=sub_ecmwf.ecmwf_ens_area(area['col2'][i],index1,index2,hs2,time_len)
    wd_mi,wd_ma,wd_75,wd_25,wd_50,wd_median=sub_ecmwf.ecmwf_ens_area(area['col2'][i],index1,index2,wd2,time_len)
    wp_mi,wp_ma,wp_75,wp_25,wp_50,wp_median=sub_ecmwf.ecmwf_ens_area(area['col2'][i],index1,index2,wp2,time_len)
    wind_vel_mi,wind_vel_ma,wind_vel_75,wind_vel_25,wind_vel_50,wind_vel_median=sub_ecmwf.ecmwf_ens_area(area['col2'][i],index1,index2,wind_vel2,time_len)
    wind_dir_mi,wind_dir_ma,wind_dir_75,wind_dir_25,wind_dir_50,wind_dir_median=sub_ecmwf.ecmwf_ens_area(area['col2'][i],index1,index2,wind_dir2,time_len)
    
    out_file=dasv+'/AREA_ecmwf_ense_'+str(area['col2'][i])+'_'+srtdate+'.dat'
    f=open(out_file,'w')
    f.write('                                                                                  HS                                                               WD                                      WP                WIND_vel\n' )
    f.write('Year Month Day Hour minimum maximum 75% 25% 50% median minimum maximum 75% 25% 50 median minimum maximum 75% 25% 50 median minimum maximum 75% 25% 50 median\n')        
    
    for j in range(len(time[:])):
        nt=Time(time[j]*3600,format='unix',scale='utc')
        nt1=st1.unix+nt.unix
        nt2=Time(nt1,format='unix',scale='utc')
        yy=nt2.iso[0:4]
        mm=nt2.iso[5:7]
        dd=nt2.iso[8:10]
        hh=nt2.iso[11:13]
        f.write('%s %s %s %s %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f %3.4f\n' \
                %(yy,mm,dd,hh,hs_mi[j],hs_ma[j],hs_75[j],hs_25[j],hs_50[j],hs_median[j],\
                  wd_mi[j],wd_ma[j],wd_75[j],wd_25[j],wd_50[j],wd_median[j],\
                  wp_mi[j],wp_ma[j],wp_75[j],wp_25[j],wp_50[j],wp_median[j],\
                  wind_vel_mi[j],wind_vel_ma[j],wind_vel_75[j],wind_vel_25[j],wind_vel_50[j],wind_vel_median[j],\
                 wind_dir_mi[j],wind_dir_ma[j],wind_dir_75[j],wind_dir_25[j],wind_dir_50[j],wind_dir_median[j]))
    f.close()    
    
    
    
    
    
