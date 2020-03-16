#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:06:18 2020

@author: Joseph Martin, University of Victoria 
with contributions from Dr. Michael Sigmond, Canadian Centre for Climate Modelling and Analysis

"""

from netCDF4 import Dataset
import math as math
import nc as nc
import numpy as np
import array as ar
import matplotlib.colors
import matplotlib.pyplot as plt
import calendar
import rms_plots as rpl
import CCMA_plot

#Select CANSIPS v1 or v2

CANSIPS="v2"

#Select OLD or NEW
versions=["OLD","NEW"]

#Set plot style
pstyle = "pc"
if pstyle != "cf" and pstyle != "pc":
    print("Select a plot style (cf or pc)")
    
for version in versions:    
    #Select SIE or SIA
    metric="SIA"
    
    obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', str.lower(metric)).squeeze()
    
    #Select Year Range
    years=range(1980,2010)
    
    sim=np.zeros((12,12,len(years)))
    
    #Shape observations
    obs=np.delete((np.reshape(obsin, (32,12)).transpose()),range(min(years)-1979),1)
    
    
    
    
    #Build model array
    for months in range(1,13):
    	if months<10:
    		m="0%s" % months
    	else:
    		m=str(months)
    	for year in years:
                y=str(year)
                if CANSIPS=="v1":
                    var3=nc.getvar(('/home/josmarti/Data/%s/%s/%s_monthly_CanCM3_i%s%s.nc' % (version,metric,metric,y,m)),'sic').squeeze()
                if CANSIPS=="v2":
                    var3=nc.getvar(('/home/josmarti/Data/GEM-NEMO/%s/%s_monthly_GEM-NEMO_i%s%s.nc' % (metric,metric,y,m)),'sic').squeeze()
                var4=nc.getvar(('/home/josmarti/Data/%s/%s/%s_monthly_CanCM4_i%s%s.nc' % (version,metric,metric,y,m)),'sic').squeeze()
                var=np.concatenate((var3,var4), axis=1)
                extent=var*2*math.pi*6.371**2 #multiply constant to convert fraction to SIE
                #avex=np.mean(extent, axis=1)
                sim[months-1,:,year-min(years)]=np.mean(extent, axis=1) #Average across ensembles and insert into matrix
    
    
    for i in range(12): #line up all target months
        sim[i,:,:]=np.roll(sim[i,:,:],i,axis=0)
    
    
    #Calculate ACC
    ACC=np.zeros((12,12), dtype=np.ndarray)
    for init in range(12):
    	for target in range(12):
                if init<=target: 
                    ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,0:-1])[1,0]
                else: #rolls observations to realign years
                    ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,1::])[1,0]
    
    ACC = np.vstack(ACC[:, :]).astype(np.float) #ACC is an object and pcolor needs floats
    
    
    
    ACC2=np.zeros((12,12))

    for init in range(12):
    
        
        for target in range(12):
        
        
            if target>=init: #same year
            
            
                ilag=target-init
            
            
            else: #next year 
            
            
                ilag=target-init+12


            ACC2[ilag,target]=ACC[init,target] 
            
    if version=="OLD":
        old_ACC=ACC2
    
    if version=="NEW":
        new_ACC=ACC2
     
    fig, ax = plt.subplots()
    if pstyle == "pc":
        pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
        pc=rpl.add_pc(ax,range(13),range(13),ACC2,**pcparams)
    elif pstyle == "cf":
        pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar', latlon=False)
        pc=rpl.add_cf(ax,range(12),range(12),ACC2,**pcparams)
    plt.colorbar(pc)
    plt.title("%s %s ACC of forecasts from %i to %i" % (str.capitalize(version), metric, min(years),(max(years)+1)))
    ax.invert_xaxis
    ax.invert_yaxis
    ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
    ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
    ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
    ax.set_yticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])    
    plt.xlabel("Predicted Month")    
    plt.ylabel("Lead Time")      
    plt.show()
    
    
    """
    k=pd.DataFrame(ACC2,index=calendar.month_name[1:13], columns=range(1,13))
    fig, ax = plt.subplots()
    d=ax.pcolor(k)
    plt.colorbar(d)
    plt.title("%s %s ACC of forecasts from %i to %i" % (str.capitalize(version), metric, min(years),(max(years)+1)))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    #ax.xaxis.set_minor_locator(ticker.FixedLocator((range(1,13)+0.5)))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(str(range(1,13))))
    ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)
    ax.set_yticklabels(k.index)
    plt.xlabel("Lead Time")
    plt.ylabel("Initialization month")
    plt.show() """
    
ACC=new_ACC-old_ACC

fig, ax = plt.subplots()
if pstyle == "pc":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
    pc=rpl.add_pc(ax,range(13),range(13),ACC,**pcparams)
elif pstyle == "cf":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar',latlon=False)
    pc=rpl.add_cf(ax,range(1,13),range(1,13),ACC,**pcparams)
plt.colorbar(pc)
plt.title("Difference in ACC of forecasts from %i to %i" % (min(years),(max(years)+1)))
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
ax.set_yticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])    
plt.xlabel("Predicted Month")    
plt.ylabel("Lead Time")      
plt.show()