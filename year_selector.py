#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 08:29:19 2020

@author: Joseph Martin, University of Victoria
"""

import nc as nc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



years=range(1990,2300)

siv=pd.DataFrame(0, index=range(12), columns=years)
gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')

for year in years:
    var=nc.getvar(('/home/josmarti/Data/CanCM4_control/sc_dhfp1e_e001_%s_m01_%s_m12_sicn.nc' % (year,year)),'sicn').squeeze()
    sit=var/913.0 #SIT=SIC/913 kg/m3
    siv_geo_nh=np.delete(np.multiply(sit,gridpoint), range(32), axis=1)
    siv[year]=np.sum(np.sum(siv_geo_nh, axis=2), axis=1, keepdims=True)
    
plt.figure(figsize=[12,9])
plt.plot(years, siv.mean(axis=0))
plt.title('Annual Mean Sea Ice Volume of Control Run \n')
plt.xlabel('Years')
plt.ylabel('SIV')
    
siv_mean=siv.mean(axis=1)

siv_anoms=siv.subtract(siv_mean, axis=0)

siv_anoms_mean=siv_anoms.mean(axis=0).sort_values()

sextiles=np.array_split(siv_anoms_mean, 6)

selected=pd.DataFrame(index=range(len(sextiles)), columns=['Year', 'Mean of Monthly Anomalies'])

separation_margin=20
tries=range(1000)

for t in tries:
    if t%100==0:
        print ('%s percent of tries complete \n' % str((t*100/len(tries))))
    if t==max(tries):
        print ('No successful selections')
    separation=[]
    for s in range(len(sextiles)):
        pick=np.random.randint(0, high=len(sextiles[s])-1)
        selected.iloc[s,0]=sextiles[s].index[pick]
        if selected.iloc[s,0] in separation:
            break
        else:
            separation.extend(range(selected.iloc[s,0]-separation_margin,selected.iloc[s,0]+separation_margin))
            #print separation
        selected.iloc[s,1]=sextiles[s].iloc[pick]
    if s==5 and not pd.isna(selected.iloc[5,1]):
        break
print selected
print ('\nSelected with a separation of %s years after %s tries.' % (str(separation_margin),str(t)))

plt.figure(figsize=[12,9])
plt.plot(years, siv_anoms.mean(axis=0))
plt.scatter(selected['Year'], selected['Mean of Monthly Anomalies'], color='r', s=50)
plt.title('Annual Sea Ice Volume Anomalies of Control Run \n')
plt.xlabel('Years')
plt.ylabel('Anomaly')