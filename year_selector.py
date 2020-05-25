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

sic=pd.DataFrame(0, index=range(12), columns=years)

for year in years:
    var=nc.getvar(('/home/josmarti/Data/CanCM4_control/sc_dhfp1e_e001_%s_m01_%s_m12_sic.nc' % (year,year)),'sic').squeeze()
    sic[year]=np.mean(np.mean(var, axis=2), axis=1, keepdims=True)
    
sic_mean=sic.mean(axis=1)

sic_anoms=sic.subtract(sic_mean, axis=0)

sic_anoms_mean=sic_anoms.mean(axis=0).sort_values()

sextiles=np.array_split(sic_anoms_mean, 6)

selected=pd.DataFrame(index=range(len(sextiles)), columns=['Year', 'Anomaly'])

separation_margin=20
tries=range(1000)

for t in tries:
    if t%100==0:
        print ('%s percent of tries complete \n' % str(t))
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

plt.plot(years, sic_anoms.mean(axis=0))
plt.scatter(selected['Year'], selected['Anomaly'], color='r', s=50)
plt.title('Mean Monthly Sea Ice Concentration Anomalies of Control Run')
plt.xlabel('Years')
plt.ylabel('Anomaly')