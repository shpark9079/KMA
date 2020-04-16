#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 13:49:43 2019

@author: shpark
"""

from scipy.io  import netcdf
from astropy.io import ascii
import sys
import numpy as np
from astropy.time import Time
import numpy.ma as ma
import static_test

def wind_dir_statistic(index,index1,index2,uv,u,v,*param):
    
    
    min1=[]
    max1=[]
    median1=[]
    avg1=[]    
    num1=[]
    ind=[601,602,603,604,605]
    ind_du=np.where(index==ind)
    
    if (len(ind_du[0]) !=0):
        ind1=index2
    else:
        ind1=index1
    
    dumy=np.where(ind1==index)
    for i in range(u.shape[0]):
        u1=u[i,dumy[0],dumy[1]]
        v1=v[i,dumy[0],dumy[1]]
        uv1=uv[i,dumy[0],dumy[1]]
        u2=u1[~u1.mask]        
        v2=v1[~v1.mask]
        uv2=uv1[~uv1.mask]
        dumy_min=np.where(uv2==np.min(uv2))
        dumy_max=np.where(uv2==np.max(uv2))
        if ((len(dumy_min[0])>0) | (len(dumy_max[0]>0))):
            u_mi=u2[dumy_min]
            v_mi=v2[dumy_min]
            wd_mi=np.arctan2(v_mi,u_mi)*(180/np.pi)
            wd_mi1=(450-wd_mi)%360
            u_ma=u2[dumy_max]
            v_ma=v2[dumy_max]
            wd_ma=np.arctan2(v_ma,u_ma)*(180/np.pi)
            wd_ma1=(450-wd_ma)%360
            u_me=np.median(u2)
            v_me=np.median(v2)
            u_av=np.mean(u2)
            v_av=np.mean(v2)
            wd_me=np.arctan2(v_me,u_me)*(180/np.pi)
            wd_me1=(450-wd_me)%360
            wd_av=np.arctan2(v_av,u_av)*(180/np.pi)
            wd_av1=(450-wd_av)%360
        else:
            wd_mi1=np.nan
            wd_ma1=np.nan
            wd_me1=np.nan
            wd_av1=np.nan
                    
        num2=ma.count_masked(uv1)
        num3=len(uv1)-num2
        num1.append(num3)
        min1.append(wd_mi1[0])
        max1.append(wd_ma1[0])
        median1.append(wd_me1)
        avg1.append(wd_av1)        
        
    return num1,min1,max1,median1,avg1        
#args=sys.argv[1:]
#daio=args[0]
#dain=args[1]
#dasv=args[2]
#srtdate=args[3]
#netcdf_file=args[4]

dain='/data1/2019/KMA/KMA_2019/SFP_1/1.1/model/OAWN/DAIN'
daio='/data1/2019/KMA/KMA_2019/SFP_1/1.1/model/OAWN/DAIO'
dasv='/data1/2019/KMA/KMA_2019/SFP_1/1.1/model/OAWN/DASV'
srtdate='2019050700'
netcdf_file=daio+'/rww3.'+srtdate+'.nc'
#base_dir='/data1/2019/KMA/KMA_2019/SFP_1/1.1/model'
#area_index=base_dir+'/area_index.dat'
area_index=dain+'/area_index1.dat'

ecmwf_nc=netcdf.netcdf_file(daio+'/ecmwf_ense1.nc','r')
#ecmwf_nc=netcdf.netcdf_file(daio+'/'+netcdf_file,'r')
hs=ecmwf_nc.variables['swh']
wd=ecmwf_nc.variables['mwd']
wp=ecmwf_nc.variables['mwp']
wind_vel=ecmwf_nc.variables['wind']
wind_dir=ecmwf_nc.variables['dwi']
fill_v=hs._FillValue
hs1=ma.masked_where(hs[:,:,:]==hs._FillValue,hs[:,:,:])
wd1=ma.masked_where(wd[:,:,:]==wd._FillValue,wd[:,:,:])
wp1=ma.masked_where(wp[:,:,:]==wp._FillValue,wp[:,:,:])
wind_vel1=ma.masked_where(wind_vel[:,:,:]==wind_vel._FillValue,wind_vel[:,:,:])
wind_dir1=ma.masked_where(wind_dir[:,:,:]==wind_dir._FillValue,wind_dir[:,:,:])
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
ncid1=netcdf.netcdf_file(dain+'/ecmwf_only_index1.nc')
ncid2=netcdf.netcdf_file(dain+'/ecmwf_only_index2.nc')

index1=ncid1.variables['index'][:,:]
index2=ncid2.variables['index'][:,:]

for i in range(len(area['col1'])):
    hs_num1,hs_min1,hs_max1,hs_median1,hs_avg1,hs_std1=static_test.area_cal(area['col2'][i],index1,index2,hs2) 
    wd_num1,wd_min1,wd_max1,wd_median1,wd_avg1,wd_std1=static_test.area_cal(area['col2'][i],index1,index2,wd2) 
    wp_num1,wp_min1,wp_max1,wp_median1,wp_avg1,wp_std1=static_test.area_cal(area['col2'][i],index1,index2,wp) 
    vel_num1,vel_min1,vel_max1,vel_median1,vel_avg1,vel_std1=static_test.area_cal(area['col2'][i],index1,index2,wind_vel2)    
    wind_num1,wind_min1,wind_max1,wind_median1,wind_avg1=static_test.area_cal(area['col2'][i],index1,index2,wind_dir2)    
    
    #out_file=base_dir+'/out1/AREA_rww3_'+str(area['col2'][i])+'_2019050700.dat'
    out_file=dasv+'/AREA_ecmwf_only_'+str(area['col2'][i])+'_'+srtdate+'.dat'
    f=open(out_file,'w')
    f.write('YEAR Month Day Hour Hs WD WP WIND_vel WIND_dir\n')
    for j in range(len(time[:])):
        nt=Time(time[j],format='unix',scale='utc')
        yy=nt.iso[0:4]
        mm=nt.iso[5:7]
        dd=nt.iso[8:10]
        hh=nt.iso[11:13]
        f.write('%s %s %s %s %5d %3.4f %3.4f %3.4f %3.4f %5d %3.4f %3.4f %3.4f %3.4f %5d %3.4f %3.4f %3.4f %3.4f  %5d %3.4f %3.4f %3.4f %3.4f %5d %3.4f %3.4f %3.4f %3.4f\n' \
                %(yy,mm,dd,hh,hs_num1[j],hs_min1[j],hs_max1[j],hs_median1[j],hs_avg1[j],\
                  wd_num1[j],wd_min1[j],wd_max1[j],wd_median1[j],wd_avg1[j],\
                  wp_num1[j],wp_min1[j],wp_max1[j],wp_median1[j],wp_avg1[j],\
                  vel_num1[j],vel_min1[j],vel_max1[j],vel_median1[j],vel_avg1[j],\
                 wind_num1[j],wind_min1[j],wind_max1[j],wind_median1[j],wind_avg1[j]))
    f.close()       
