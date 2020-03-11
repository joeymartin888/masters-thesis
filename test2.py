#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 22:43:59 2020

@author: acrnrms
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:06:18 2020
@author: Joseph Martin, University of Victoria (edited by Michael Sigmond)
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

basedir='/home/josmarti/Data/SIE'
obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sie').squeeze()

sim=np.zeros((12,12,31))

#Shape observations
obs=np.reshape(obsin, (32,12)).transpose()




#Build model array
for months in range(1,13):
	if months<10:
		m="0%s" % months
	else:
		m=str(months)
	for years in range (1979,2009+1):
            y=str(years)
            var3=nc.getvar(('{}/SIE_monthly_CanCM3_i%s%s.nc' % (y,m)).format(basedir),'sic').squeeze()
            var4=nc.getvar(('{}/SIE_monthly_CanCM4_i%s%s.nc' % (y,m)).format(basedir),'sic').squeeze()
            var=np.concatenate((var3,var4), axis=1)
            extent=var*2*math.pi*6.371**2 #multiply constant to convert fraction to SIE
            #avex=np.mean(extent, axis=1)
            sim[months-1,:,years-1979]=np.mean(extent, axis=1) #Average across ensembles and insert into matrix


for i in range(12): #line up all target months
    sim[i,:,:]=np.roll(sim[i,:,:],i,axis=0)


#Calculate ACC
ACC=np.zeros((12,12))
for init in range(12):
	for target in range(12):
            if target>=init: #same year
              ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,0:-1])[0,1]
            else: #next year 
              ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,1::])[0,1]
            

#ACC = np.vstack(ACC[:, :]).astype(np.float) #ACC is an object and pcolor needs floats

##################Reproduction of Fig1 S13#############################3
years=np.arange(1979,2010)

fig1, axs = plt.subplots(1,1, figsize=(8,8)); fig1.subplots_adjust(right=0.6,bottom=0.5)

#yaxis
axs.set_ylabel('SIE (10$^6$ km$^2$)')
#axs.set_ylim([-2.6, 1.5])
#NSIDC
kwargs={'linewidth': 1, 'color': 'gray','label':'NSIDC'}
axs.plot(years,obs[8,0:-1].transpose(),**kwargs)
#Model lag=0
kwargs={'linewidth': 1, 'color': 'blue','label':'lag=0'}
axs.plot(years,sim[8,8,:].transpose(),**kwargs)
print ACC[8,8]
#Model lag=1
kwargs={'linewidth': 1, 'color': 'green','label':'lag=1'}
axs.plot(years,sim[7,8,:].transpose(),**kwargs)
print ACC[7,8]
#Model lag=4
kwargs={'linewidth': 1, 'color': 'red','label':'lag=3'}
axs.plot(years,sim[5,8,:].transpose(),**kwargs)
print ACC[5,8]
#Model lag=11
kwargs={'linewidth': 1, 'color': 'orange','label':'lag=11'}
axs.plot(years,sim[9,8,:].transpose(),**kwargs)
print ACC[9,8]

##################Testplot September#############################3


#Plot ACC
k=pd.DataFrame(ACC,index=calendar.month_name[1:13], columns=calendar.month_name[1:13])
fig, ax = plt.subplots()
d=ax.pcolor(k)
plt.colorbar(d)
plt.title("ACC of forecasts from 1979 to 2009")
#ax.invert_xaxis
#ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_xticklabels(k.index, rotation=90)
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)
ax.set_yticklabels(k.index)
plt.xlabel("Target month")
plt.ylabel("Initialization month")
plt.show()