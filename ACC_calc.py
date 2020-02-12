#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:06:18 2020

@author: Joseph Martin, University of Victoria
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



obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sie')

sim=np.zeros((12,12,32))

#Shape observations
obs=obsin[:,0,0,0]
obs=np.reshape(obs, (32,12)).transpose()
obs2=np.delete(obs,0,1)



#Build model array
for months in range(1,13):
	if months<10:
		m="0%s" % months
	else:
		m=str(months)
	for years in range (1980,2011):
            y=str(years)
            var3=nc.getvar(('/home/josmarti/Data/SIE/SIE_monthly_CanCM3_i%s%s.nc' % (y,m)),'sic')
            var3=var3[:,:,0,0] #strips extra dimensions (check which ones)
            var4=nc.getvar(('/home/josmarti/Data/SIE/SIE_monthly_CanCM4_i%s%s.nc' % (y,m)),'sic')
            var4=var4[:,:,0,0]
            var=np.concatenate((var3,var4), axis=1)
            extent=var*2*math.pi*6.371**2 #multiply constant to convert fraction to SIE
            #avex=np.mean(extent, axis=1)
            sim[months-1,:,years-1979]=np.mean(extent, axis=1) #Average across ensembles and insert into matrix


for i in range(12): #line up all target months
    sim[i,:,:]=np.roll(sim[i,:,:],i,axis=0)

sim2=np.delete(sim,0,2)

#Calculate ACC
ACC=np.zeros((12,12), dtype=np.ndarray)
for init in range(12):
	for target in range(12):
            ACC[init,target]=np.corrcoef(sim2[init,target,:],obs2[target,:])[1,0]
            

ACC = np.vstack(ACC[:, :]).astype(np.float) #ACC is an object and pcolor needs floats

print(ACC)
#Plot ACC
k=pd.DataFrame(ACC,index=calendar.month_name[1:13], columns=calendar.month_name[1:13])
fig, ax = plt.subplots()
ax.pcolor(k)
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_xticklabels(k.index, rotation=90)
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_yticklabels(k.index)
plt.colorbar(ACC)
plt.show()

