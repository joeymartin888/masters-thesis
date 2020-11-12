#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

@author: josmarti
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import cartopy.crs as ccrs
import rms_plots as rpl
import nc as nc
import numpy as np
import re
import calendar
import rms_utils_boot as bt
import CCMA_plot
from scipy.io import loadmat
import math as math
#import cdo; c=cdo.Cdo()

#%%

pstyle="pc"

metric='SIE'

syear=1979; eyear=2016; 

obsin=np.zeros((12,40)) 

if metric=='SIE':
    for m in ["%02d" % i for i in range(1,13)]:
            obsin[int(m)-1,:]=np.loadtxt('/home/josmarti/Data/Observations/NSIDC/north/N_%s_extent_v3.0.csv' % m,usecols=4, skiprows=1, delimiter=',', max_rows=40)
    obsin[obsin == -9999] = np.nan
if metric=='SIA':
    for m in ["%02d" % i for i in range(1,13)]:
            obsin[int(m)-1,:]=np.loadtxt('/home/josmarti/Data/Observations/NSIDC/north/N_%s_extent_v3.0.csv' % m,usecols=5, skiprows=1, delimiter=',', max_rows=40)
    obsin[obsin == -9999] = np.nan
mean_obs=np.nanmean(obsin[:,(syear-1979):(eyear-1979)], axis=1)*1e12
annual_mean_obs=np.nanmean(obsin[:,(syear-1979):(eyear-1979)], axis=0)*1e12
#gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')

PIOMAS=np.loadtxt('/home/josmarti/Data/Observations/PIOMAS_SIV_monthly_1979-2020.txt',usecols=range(1,13))
if metric=='SIV':
    mean_obs=np.mean(PIOMAS[(syear-1979):(eyear-1979)], axis=0)*1e12

years=[2312, 2435, 2413, 2356, 2387, 2332]
ensembles=["%03d" % i for i in range(1,13)]
#ensembles=["001"]
run_length=1 #years
daily_matrix=np.zeros((6,12,1095))
days=365*3
smonth=["%02d" % i for i in range(1,12,2)]
emonth=["%02d" % i for i in range(0,12,2)]
emonth[0]="12"

monthly_matrix=np.zeros((6,72,36))

for i in range(6):
    if i == 0:
        jd=0
    else:
        jd=1
    for e in ensembles:
        for year in years:
            if metric=='SIE':
                obs=annual_mean_obs
                geosum1=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PP_run_full_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+2+jd),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
                #geosum2=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PPrun_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+1),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
                #geosum3=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PPrun_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+2),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
            if metric=='SIA':
                obs=annual_mean_obs
                var1=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year),str(year))), 'sicn'), gridpoint)
                var2=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year+1),str(year+1))), 'sicn'), gridpoint)
                var3=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year+2),str(year+2))), 'sicn'), gridpoint)
            if metric=='SIV':
                obs=np.mean(PIOMAS[(syear-1979):(eyear-1979)], axis=1)*1e12
                var1=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year),str(year))), 'sic'), gridpoint/913)
                var2=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+1),str(year+1))), 'sic'), gridpoint/913)
                var3=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+2),str(year+2))), 'sic'), gridpoint/913)
            #var=nc.getvar('/home/josmarti/Data/PPrun_Jan/sc_dhfp1e_e%s_i%s_m01_daily_sicn.nc4' % (e,str(year)),'sicn')
            if metric != 'SIE':
                var1nh=np.delete(var1, range(32), axis=1)
                var2nh=np.delete(var2, range(32), axis=1)
                var3nh=np.delete(var3, range(32), axis=1)
                geosum1=np.sum(np.sum(var1nh, axis=1), axis=1)
                geosum2=np.sum(np.sum(var2nh, axis=1), axis=1)
                geosum3=np.sum(np.sum(var3nh, axis=1), axis=1)
            #daily_matrix[years.index(year),ensembles.index(e),:]=np.sum(np.sum(var[:,32:64,:], axis=1), axis=1)
            if year==2413 and e=='011' and i==0: #i241501 e11 seems to have an extra timestep
                geosum1=geosum1[0:-1]
            monthly_matrix[i,(12*years.index(year)+ensembles.index(e)),:]=geosum1

#%%
for j in range(6):            
    monthly_matrix[j,:,:]=np.reshape(np.transpose(loadmat('/home/josmarti/Downloads/SIE_16.mat')['metric_ensemble'][:,j,:,:], (0,2,1)), (72,36))

#%%
model=np.zeros((6,72,36))
truth=np.zeros((6,72,36))

ACC=np.zeros((12,36))
boot_temp=np.zeros((12,36), dtype=np.ndarray)
boot=np.zeros((12,36), dtype=np.ndarray)

for y in range(6):
    for e in range(12):
        model[:,(y*12+e),:]=np.mean(np.delete(monthly_matrix[:,y*12:y*12+12,:], e, axis=1), axis=1)
        #truth[:,e,:]=monthly_matrix[]

truth=monthly_matrix

for i in range(0,12,2):
    for m in range(36):
        #model=np.repeat(np.mean(monthly_matrix[:,:,m], axis=1),12)
        #truth=np.reshape(monthly_matrix[:,:,m],72)
        if m<12:
            ACC[m-i,m]=np.corrcoef(model[(i/2),:,m],truth[(i/2),:,m])[1,0]
            boot_temp[m-i,m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[m-i,m]=bt.calc_boot_stats(boot_temp[m-i,m],sides=1,pval_threshold=0.05)[1]
        elif m<24:
            ACC[m-12-i,m]=np.corrcoef(model[(i/2),:,m],truth[(i/2),:,m])[1,0]
            boot_temp[m-12-i,m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[m-12-i,m]=bt.calc_boot_stats(boot_temp[m-12-i,m],sides=1,pval_threshold=0.05)[1]
        else:
            ACC[m-24-i,m]=np.corrcoef(model[(i/2),:,m],truth[(i/2),:,m])[1,0]
            boot_temp[m-24-i,m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[m-24-i,m]=bt.calc_boot_stats(boot_temp[m-24-i,m],sides=1,pval_threshold=0.05)[1]

#%%         

#ACC=data['ACC_year_target']
jet3=colors.ListedColormap(loadmat('/home/josmarti/Downloads/cmap_jet3.mat')['cmap'], name='jet3') #Mitch's colors


fig, ax = plt.subplots(figsize=(4,9))
if pstyle == "pc":
    #pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar') #CCCMA
    pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)  #Mitch's
    pc=rpl.add_pc(ax,range(13),range(37),ACC.transpose(),**pcparams)
elif pstyle == "cf":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar',latlon=False)
    pc=rpl.add_cf(ax,range(1,13),range(1,37),ACC,**pcparams)
for month in range(len(boot[:,0])):
            for lead in range(len(boot[0,:])):
                if boot[month,lead]>0:   #MUST BE CHANGED BACK to >=!!!!!!
                    plt.scatter((month+0.5),(lead+0.5), color = 'black', s=20)
plt.title("Forecast skill for %s of PP_run" % metric)
plt.colorbar(pc)
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
ax.set_yticks(np.arange(36)+0.5)    
ax.set_yticklabels(range(36))  
plt.xlabel("Predicted Month")    
plt.ylabel("Lead (Months)")      
plt.show()            
        