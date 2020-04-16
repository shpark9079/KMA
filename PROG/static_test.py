#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:30:52 2019

@author: shpark
"""

import numpy as np
import wind_dir_cal
import numpy.ma as ma
from astropy.time import Time

def area_cal(index,index1,index2,vari,*param):
    
    min1=[]
    max1=[]
    median1=[]
    avg1=[]
    std1=[]
    num1=[]
    ind=[601,602,603,604,605]
    ind_du=np.where(index==ind)
    
    if (len(ind_du[0]) !=0):
        ind1=index2
    else:
        ind1=index1
    
    dumy=np.where(ind1==index)
    for i in range(vari.shape[0]):
          min2=ma.min(vari[i,dumy[0],dumy[1]])
          max2=ma.max(vari[i,dumy[0],dumy[1]])
          median2=ma.median(vari[i,dumy[0],dumy[1]])
          avg2=ma.mean(vari[i,dumy[0],dumy[1]])
          std2=ma.std(vari[i,dumy[0],dumy[1]])
          num2=ma.count_masked(vari[i,dumy[0],dumy[1]])
          num3=len(vari[i,dumy[0],dumy[1]])-num2
          num1.append(num3)
          min1.append(min2)
          max1.append(max2)
          median1.append(median2)
          avg1.append(avg2)
          std1.append(std2)
    return num1,min1,max1,median1,avg1,std1

def area_rdwa_cal(index,index1,index2,vari,*param):
    
    
    ind=[601,602,603,604,605]
    ind_du=np.where(index==ind)
    
    if (len(ind_du[0]) !=0):
        ind1=index2
    else:
        ind1=index1
    
    dumy=np.where(ind1==index)
    min2=ma.min(vari[dumy[0],dumy[1]])
    max2=ma.max(vari[dumy[0],dumy[1]])
    median2=ma.median(vari[dumy[0],dumy[1]])
    avg2=ma.mean(vari[dumy[0],dumy[1]])
    std2=ma.std(vari[dumy[0],dumy[1]])
    num2=ma.count_masked(vari[dumy[0],dumy[1]])
    num3=len(vari[dumy[0],dumy[1]])-num2
    
    return num3,min2,max2,median2,avg2

def area_cal1(index,index1,vari,*param):
    
    min1=[]
    max1=[]
    median1=[]
    avg1=[]
    std1=[]
    num1=[]
    ind1=index1
        
    dumy=np.where(ind1==index)
    for i in range(vari.shape[0]):
          min2=ma.min(vari[i,dumy[0],dumy[1]])
          max2=ma.max(vari[i,dumy[0],dumy[1]])
          median2=ma.median(vari[i,dumy[0],dumy[1]])
          avg2=ma.mean(vari[i,dumy[0],dumy[1]])
          std2=ma.std(vari[i,dumy[0],dumy[1]])
          num2=ma.count_masked(vari[i,dumy[0],dumy[1]])
          num3=len(vari[i,dumy[0],dumy[1]])-num2
          num1.append(num3)
          min1.append(min2)
          max1.append(max2)
          median1.append(median2)
          avg1.append(avg2)
          std1.append(std2)
    return num1,min1,max1,median1,avg1

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

def wind_dir_statistic1(index,index1,u,v,*param):
    
    min1=[]
    max1=[]
    median1=[]
    avg1=[]
    std1=[]
    num1=[]
    ind1=index1
    
    dumy=np.where(ind1==index)
    for i in range(u.shape[0]):
        u1=u[i,dumy[0],dumy[1]]
        v1=v[i,dumy[0],dumy[1]]
        u2=u1[~u1.mask]        
        v2=v1[~v1.mask]
        wind_dir=wind_dir_cal.wind_dir_cal(u2,v2)        
        min2=ma.min(wind_dir)
        max2=ma.max(wind_dir)
        median2=ma.median(wind_dir)
        avg2=ma.mean(wind_dir)
        std2=ma.std(wind_dir)
        num2=ma.count_masked(wind_dir)
        num3=len(wind_dir)-num2
        num1.append(num3)
        min1.append(min2)
        max1.append(max2)
        median1.append(median2)
        avg1.append(avg2)
        std1.append(std2)
        
    return num1,min1,max1,median1,avg1,std1        
        
def rdwa_marinezone(dasv,index,lon,lat,var,date,*param):
    
    loc_lon=index['col3']
    loc_lat=index['col2']
    xx1,yy1=np.meshgrid(lon,lat)
    outfile=dasv+'/marinezone_rdwa_'+date+'.dat'
    f=open(outfile,'w')
    f.write('zone_num minimum maximum median average std\n')
    for i in range(len(index['col1'])):
        
        dumy=np.where((loc_lon[i]<=xx1)&(loc_lon[i]+0.5>=xx1)&(loc_lat[i]<=yy1)&(loc_lat[i]+0.5>=yy1))
        var1=var[dumy[0],dumy[1]]
        num=np.isnan(var1).sum()
        num1=len(var1)-num
        var_min=np.nanmin(var1)
        var_max=np.nanmax(var1)
        var_me=np.nanmedian(var1)
        var_av=np.nanmean(var1)
        var_st=np.nanstd(var1)
        f.write('%5d %6d %3.4f %3.4f %3.4f %3.4f %3.4f\n' %(index['col1'][i],num1,var_min,var_max,var_me,var_av,var_st))
    f.close()
        
        
def rww3_marinezone(dasv,index,lon,lat,hs,wd,wp,u,v,time,*param):
    
    loc_lon=index['col3']
    loc_lat=index['col2']
    xx1,yy1=np.meshgrid(lon,lat)
    
    for t in range(len(time)):
        nt=Time(time[t],format='unix',scale='utc')
        yy=nt.iso[0:4]
        mm=nt.iso[5:7]
        dd=nt.iso[8:10]
        hh=nt.iso[11:13]
        ndate=str(yy)+str(mm)+str(dd)+str(hh)
        outfile=dasv+'/marinezone_rww3_'+ndate+'.dat'
        f=open(outfile,'w')
        f.write('zone_num minimum maximum median average std\n')
        hs_var=hs[t,:,:]
        wd_var=wd[t,:,:]
        wp_var=wp[t,:,:]
        u_var=u[t,:,:]
        v_var=v[t,:,:]
        for i in range(len(index['col1'])):
            
            dumy=np.where((loc_lon[i]<=xx1)&(loc_lon[i]+0.5>=xx1)&(loc_lat[i]<=yy1)&(loc_lat[i]+0.5>=yy1))
            hs_var1=hs_var[dumy[0],dumy[1]]
            wd_var1=wd_var[dumy[0],dumy[1]]
            wp_var1=wp_var[dumy[0],dumy[1]]
            u_var1=u_var[dumy[0],dumy[1]]
            v_var1=v_var[dumy[0],dumy[1]]
            hs1=hs_var1[~hs_var1.mask]
            wd1=wd_var1[~wd_var1.mask]
            wp1=wp_var1[~wp_var1.mask]
            u1=u_var1[~u_var1.mask]
            v1=v_var1[~v_var1.mask]
            if ((len(u1) == 0) & (len(v1)==0) & (len(hs1)==0)):
                hs_num1=0
                wd_num1=0
                wp_num1=0
                u_num1=0
                v_num1=0
                hs_min=np.nan
                hs_max=np.nan
                hs_me=np.nan
                hs_av=np.nan
                hs_st=np.nan
                wd_min=np.nan
                wd_max=np.nan
                wd_me=np.nan
                wd_av=np.nan
                wd_st=np.nan
                wp_min=np.nan
                wp_max=np.nan
                wp_me=np.nan
                wp_av=np.nan
                wp_st=np.nan
                u_min=np.nan
                u_max=np.nan
                u_me=np.nan
                u_av=np.nan
                u_st=np.nan
                v_min=np.nan
                v_max=np.nan
                v_me=np.nan
                v_av=np.nan
                v_st=np.nan
                f.write('%5d %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s\n' \
                        %(index['col1'][i],hs_num1,hs_min,hs_max,hs_me,hs_av,hs_st,\
                          wd_num1,wd_min,wd_max,wd_me,wd_av,wd_st,
                          wp_num1,wp_min,wp_max,wp_me,wp_av,wp_st,
                          u_num1,u_min,u_max,u_me,u_av,u_st,v_num1,v_min,v_max,v_me,v_av,v_st))
            else:
                uv=np.sqrt(u1**2+v1**2)
                uv_min=np.nanmin(uv)
                uv_max=np.nanmax(uv)
                dumy_min=np.where(uv==uv_min)
                dumy_max=np.where(uv==uv_max)           
                
            
                hs_num=ma.count_masked(hs_var1)
                hs_num1=len(hs_var1)-hs_num
                wd_num=ma.count_masked(wd_var1)
                wd_num1=len(wd_var1)-wd_num
                wp_num=ma.count_masked(wp_var1)
                wp_num1=len(wp_var1)-wp_num
                u_num=ma.count_masked(u_var1)
                u_num1=len(u_var1)-u_num
                v_num=ma.count_masked(v_var1)
                v_num1=len(v_var1)-v_num                      
                hs_min=ma.min(hs1)
                hs_max=ma.max(hs1)
                hs_me=ma.median(hs1)
                hs_av=ma.mean(hs1)
                hs_st=ma.std(hs1)
                
                if ( (len(hs1) > 0 ) & (len(wd1)==0)):
                    wd_min=0
                    wd_max=0
                    wd_me=0
                    wd_av=0
                    wd_st=0
                else:
                    wd_min=ma.min(wd1)
                    wd_max=ma.max(wd1)
                    wd_me=ma.median(wd1)
                    wd_av=ma.mean(wd1)
                    wd_st=ma.std(wd1)

                wp_min=ma.min(wp1)
                wp_max=ma.max(wp1)
                wp_me=ma.median(wp1)
                wp_av=ma.mean(wp1)
                wp_st=ma.std(wp1)                    
                
                u_min=u1[dumy_min[0]]
                u_max=u1[dumy_max[0]]
                u_me=ma.median(u1)
                u_av=ma.mean(u1)
                u_st=ma.std(u1)
            
                v_min=v1[dumy_min[0]]
                v_max=v1[dumy_max[0]]
                v_me=ma.median(v1)
                v_av=ma.mean(v1)
                v_st=ma.std(v1)
                
                f.write('%5d %6d %3.4f %3.4f %3.4f %3.4f %3.4f %6d %3.4f %3.4f %3.4f %3.4f %3.4f %6d %3.4f %3.4f %3.4f %3.4f %3.4f  %6d  %3.4f %3.4f %3.4f %3.4f %3.4f  %6d %3.4f %3.4f %3.4f %3.4f %3.4f\n' \
                        %(index['col1'][i],hs_num1,hs_min,hs_max,hs_me,hs_av,hs_st,\
                          wd_num1,wd_min,wd_max,wd_me,wd_av,wd_st,\
                          wp_num1,wp_min,wp_max,wp_me,wp_av,wp_st,\
                          u_num1,u_min[0],u_max[0],u_me,u_av,u_st,v_num1,v_min[0],v_max[0],v_me,v_av,v_st))
            
        f.close()    

def ecmwf_marinezone(dasv,index,lon,lat,hs,wd,wp,u,v,time,*param):
    
    loc_lon=index['col3']
    loc_lat=index['col2']
    xx1,yy1=np.meshgrid(lon,lat)
    
    for t in range(len(time)):
        nt=Time(time[t],format='unix',scale='utc')
        yy=nt.iso[0:4]
        mm=nt.iso[5:7]
        dd=nt.iso[8:10]
        hh=nt.iso[11:13]
        ndate=str(yy)+str(mm)+str(dd)+str(hh)
        outfile=dasv+'/marinezone_ecmwf_only_'+ndate+'.dat'
        f=open(outfile,'w')
        f.write('zone_num minimum maximum median average std\n')
        hs_var=hs[t,:,:]
        wd_var=wd[t,:,:]
        wp_var=wp[t,:,:]
        u_var=u[t,:,:]
        v_var=v[t,:,:]
        for i in range(len(index['col1'])):
            
            dumy=np.where((loc_lon[i]<=xx1)&(loc_lon[i]+0.5>=xx1)&(loc_lat[i]<=yy1)&(loc_lat[i]+0.5>=yy1))
            hs_var1=hs_var[dumy[0],dumy[1]]
            wd_var1=wd_var[dumy[0],dumy[1]]
            wp_var1=wp_var[dumy[0],dumy[1]]
            u_var1=u_var[dumy[0],dumy[1]]
            v_var1=v_var[dumy[0],dumy[1]]
            hs1=hs_var1[~hs_var1.mask]
            wd1=wd_var1[~wd_var1.mask]
            wp1=wp_var1[~wp_var1.mask]
            u1=u_var1[~u_var1.mask]
            v1=v_var1[~v_var1.mask]
            if ((len(u1) == 0) & (len(v1)==0) & (len(hs1)==0)):
                hs_num1=0
                wd_num1=0
                wp_num1=0
                u_num1=0
                v_num1=0
                hs_min=np.nan
                hs_max=np.nan
                hs_me=np.nan
                hs_av=np.nan
                hs_st=np.nan
                wd_min=np.nan
                wd_max=np.nan
                wd_me=np.nan
                wd_av=np.nan
                wd_st=np.nan
                wp_min=np.nan
                wp_max=np.nan
                wp_me=np.nan
                wp_av=np.nan
                wp_st=np.nan
                u_min=np.nan
                u_max=np.nan
                u_me=np.nan
                u_av=np.nan
                u_st=np.nan
                v_min=np.nan
                v_max=np.nan
                v_me=np.nan
                v_av=np.nan
                v_st=np.nan
                f.write('%5d %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s %6d %s %s %s %s %s\n' \
                        %(index['col1'][i],hs_num1,hs_min,hs_max,hs_me,hs_av,hs_st,\
                          wd_num1,wd_min,wd_max,wd_me,wd_av,wd_st,
                          wp_num1,wp_min,wp_max,wp_me,wp_av,wp_st,
                          u_num1,u_min,u_max,u_me,u_av,u_st,v_num1,v_min,v_max,v_me,v_av,v_st))
            else:
                uv=np.sqrt(u1**2+v1**2)
                uv_min=np.nanmin(uv)
                uv_max=np.nanmax(uv)
                dumy_min=np.where(uv==uv_min)
                dumy_max=np.where(uv==uv_max)           
                
            
                hs_num=ma.count_masked(hs_var1)
                hs_num1=len(hs_var1)-hs_num
                wd_num=ma.count_masked(wd_var1)
                wd_num1=len(wd_var1)-wd_num
                wp_num=ma.count_masked(wp_var1)
                wp_num1=len(wp_var1)-wp_num
                u_num=ma.count_masked(u_var1)
                u_num1=len(u_var1)-u_num
                v_num=ma.count_masked(v_var1)
                v_num1=len(v_var1)-v_num                      
                hs_min=ma.min(hs1)
                hs_max=ma.max(hs1)
                hs_me=ma.median(hs1)
                hs_av=ma.mean(hs1)
                hs_st=ma.std(hs1)
                
                if ( (len(hs1) > 0 ) & (len(wd1)==0)):
                    wd_min=0
                    wd_max=0
                    wd_me=0
                    wd_av=0
                    wd_st=0
                else:
                    wd_min=ma.min(wd1)
                    wd_max=ma.max(wd1)
                    wd_me=ma.median(wd1)
                    wd_av=ma.mean(wd1)
                    wd_st=ma.std(wd1)

                wp_min=ma.min(wp1)
                wp_max=ma.max(wp1)
                wp_me=ma.median(wp1)
                wp_av=ma.mean(wp1)
                wp_st=ma.std(wp1)                    
                
                u_min=u1[dumy_min[0]]
                u_max=u1[dumy_max[0]]
                u_me=ma.median(u1)
                u_av=ma.mean(u1)
                u_st=ma.std(u1)
            
                v_min=v1[dumy_min[0]]
                v_max=v1[dumy_max[0]]
                v_me=ma.median(v1)
                v_av=ma.mean(v1)
                v_st=ma.std(v1)
                
                f.write('%5d %6d %3.4f %3.4f %3.4f %3.4f %3.4f %6d %3.4f %3.4f %3.4f %3.4f %3.4f %6d %3.4f %3.4f %3.4f %3.4f %3.4f  %6d  %3.4f %3.4f %3.4f %3.4f %3.4f  %6d %3.4f %3.4f %3.4f %3.4f %3.4f\n' \
                        %(index['col1'][i],hs_num1,hs_min,hs_max,hs_me,hs_av,hs_st,\
                          wd_num1,wd_min,wd_max,wd_me,wd_av,wd_st,\
                          wp_num1,wp_min,wp_max,wp_me,wp_av,wp_st,\
                          u_num1,u_min[0],u_max[0],u_me,u_av,u_st,v_num1,v_min[0],v_max[0],v_me,v_av,v_st))
            
        f.close()    

def glosea5_marinezone(dasv,index,lon,lat,var,time,nt,*param):
    
        
    loc_lon=index['col3']
    loc_lat=index['col2']
    for t in range(len(time[:])):        
        tt=Time(nt,format='iso',scale='utc')
        tt1=tt.unix+(t*10800)
        st=Time(tt1,format='unix',scale='utc')
        yy=st.iso[0:4]
        mm=st.iso[5:7]
        dd=st.iso[8:10]
        hh=st.iso[11:13]
        ndate=str(yy)+str(mm)+str(dd)+str(hh)
        outfile=dasv+'/marinezone_glosea5_'+ndate+'.dat'
        f=open(outfile,'w')
        f.write('zone_num minimum maximum median average std\n')
        for i in range(len(index['col1'])):
        
            dumy=np.where((loc_lon[i]<=lon)&(loc_lon[i]+0.5>=lon)&(loc_lat[i]<=lat)&(loc_lat[i]+0.5>=lat))
            var1=var[t,dumy[0],dumy[1]]
            num=np.isnan(var1).sum()
            num1=len(var1)-num
            var_min=np.nanmin(var1)
            var_max=np.nanmax(var1)
            var_me=np.nanmedian(var1)
            var_av=np.nanmean(var1)
            var_st=np.nanstd(var1)
            f.write('%5d %6d %3.4f %3.4f %3.4f %3.4f %3.4f\n' %(index['col1'][i],num1,var_min,var_max,var_me,var_av,var_st))
        f.close()    
