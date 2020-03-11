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
import pandas as pd
import calendar


obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sie').squeeze()

#Select OLD or NEW
version="OLD"

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
            var3=nc.getvar(('/home/josmarti/Data/%s/SIE/SIE_monthly_CanCM3_i%s%s.nc' % (version,y,m)),'sic').squeeze()
            var4=nc.getvar(('/home/josmarti/Data/%s/SIE/SIE_monthly_CanCM4_i%s%s.nc' % (version,y,m)),'sic').squeeze()
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

#print(ACC)
#Plot ACC
k=pd.DataFrame(ACC,index=calendar.month_name[1:13], columns=calendar.month_name[1:13])
fig, ax = plt.subplots()
d=ax.pcolor(k)
plt.colorbar(d)
plt.title("%s ACC of forecasts from %i to %i" % (str.capitalize(version), min(years),(max(years)+1)))
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_xticklabels(k.index, rotation=90)
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_yticklabels(k.index)
plt.xlabel("Target month")
plt.ylabel("Initialization month")
plt.show()

