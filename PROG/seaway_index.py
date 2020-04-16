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
from astropy.time import Time

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
    if (var[lat_loc[0],lon_loc[0]] <0):
        return lat_loc[0],lon_loc[0]
    else:
        for i in range(0,2):
            for j in range(0,2):
                if (var[lat_loc[0]+i,lon_loc[0]+j] < 0):
                    return lat_loc[0]+i,lon_loc[0]+j
                


#args=sys.argv[1:]
#base_dir=args[0]
#infile=args[1]
#srtdate=args[2]

base_dir='/data1/2020/KMA/SFR-07'

infile='ecmwf_only_total.nc'
#srtdate='2020021005'

ncid=netcdf.netcdf_file(base_dir+'/DAIO/'+infile,'r')

lon=ncid.variables['longitude'][:]
lat=ncid.variables['latitude'][:]
hs=ncid.variables['swh']
seaway_index=base_dir+'/DAIN/seaway_index.dat'
hs1=(hs[0,:,:]*hs.scale_factor)+hs.add_offset
area=ascii.read(seaway_index)
model_lon=lon
model_lat=lat
var=hs1
f=open(base_dir+'/DAIN/seaway_only_index.dat','w')
for k in range(len(area['col1'])):
    lon_dis=model_lon[1]-model_lon[0]
    lat_dis=abs(model_lat[1]-model_lat[0])
    lon=area['col3'][k]
    lat=area['col2'][k]
    
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
    f.write('%3d %3.4f %3.4f %4d %4d\n' % (k,lat,lon,lat_loc[0][0],lon_loc[0][0]))
    
f.close()    
#    if (var[lat_loc[0],lon_loc[0]] <0):
#        return lat_loc[0],lon_loc[0]
#    else:
#        for i in range(0,2):
#            for j in range(0,2):
#                if (var[lat_loc[0]+i,lon_loc[0]+j] < 0):
#                    return lat_loc[0]+i,lon_loc[0]+j
#time=ncid.variables['time'][:]
#hs=ncid.variables['WVHGT_surface'][:]
#wd=ncid.variables['WVDIR_surface'][:]
#wp=ncid.variables['PWPER_surface'][:]
#u=ncid.variables['UGRD_surface'][:]
#v=ncid.variables['VGRD_surface'][:]


#col_name=['stnid','lat','lon','name']
#
#yy=srtdate[0:4]
#mon=srtdate[4:6]
#dd=srtdate[6:8]
#hh=srtdate[8:10]
#
#st=yy+'-'+mon+'-'+dd+' '+hh+':00:00'
#stime=Time(st,format='iso',scale='utc')
#outfile1=base_dir+'/KWW3-loc.dat'
#f1=open(outfile1,'w')
#for loc in glob.glob(base_dir+'/DAIN/*.dat'):
##for loc in glob.glob(base_dir+'/DAIN/sea_buoy_comos.dat'):
#    loc1=pd.read_csv(loc,header=None,sep='\t')
#    loc1.columns=col_name
#    for i in range(len(loc1['stnid'])):
##    for i in range(0,12):
#        lon1=loc1['lon'][i]
#        lat1=loc1['lat'][i]
#        stn_id=str(loc1['stnid'][i])
#        name=loc1['name'][i]
#        loc_i,loc_j=fine_loc(lon,lat,lon1,lat1,hs[0,:,:])        
#        hs1=hs[:,loc_i,loc_j]            
#        wd1=wd[:,loc_i,loc_j]
#        wp1=wp[:,loc_i,loc_j]
#        u1=u[:,loc_i,loc_j]
#        v1=v[:,loc_i,loc_j]
#        outfile=base_dir+'/DASV/KWW3-'+name+'-'+stn_id+'-'+srtdate+'.dat'
#        
#        f=open(outfile,'w')
#        
#        f1.write('%s %3.4f %3.4f %s\n' % (stn_id,lat[loc_i][0],lon[loc_j][0],name))
#        f.write('model_time  fct WH WD WP U V\n')
#        for i in range(len(hs1)):
#            ft=Time(time[i],format='unix',scale='utc')
#            ftime=ft.iso[0:4]+ft.iso[5:7]+ft.iso[8:10]+ft.iso[11:13]
##            nyy=ft.iso[0:4]
##            nmon=ft.iso[5:7]
##            ndd=ft.iso[8:10]
##            nhh=ft.iso[11:14]
#            f.write('%s %s %3.4f %3.4f %3.4f %3.4f %3.4f\n' % (srtdate,ftime,hs1[i][0],wd1[i][0],wp1[i][0],u1[i][0],v1[i][0]))
#        
#        
#        f.close()
