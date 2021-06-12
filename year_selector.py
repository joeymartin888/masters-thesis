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

#regionlabs=['0land','1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15OTHER']
region="3BAR" 

if region[1].isdigit():
    r=int(region[0:2])
else:
    r=int(region[0])

mask=nc.getvar('/home/josmarti/Data/iceregions_128_64.nc', 'REG').squeeze()
mask[mask != float(r)]=0
mask[mask == float(r)]=1

years=range(2300,2449)

siv=pd.DataFrame(0, index=range(12), columns=years)
gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')

for year in years:
    var=nc.getvar(('/home/josmarti/Data/CanCM4_control/sc_dhfp1e_e001_%s_m01_%s_m12_sic.nc' % (year,year)),'sic').squeeze()
    sit=var/913.0 #SIT=SIC/913 kg/m3
    sia=nc.getvar(('/home/josmarti/Data/CanCM4_control/sc_dhfp1e_e001_%s_m01_%s_m12_sicn.nc' % (year,year)),'sicn').squeeze()
    siv_geo_nh=np.delete(np.multiply(sit,gridpoint), range(32), axis=1)
    if r!=0:
        siv_geo_nh=np.multiply(siv_geo_nh,np.delete(mask, range(32), axis=0))
    siv.iloc[0,(year-min(years))]=np.sum(siv_geo_nh)
    #siv[year]=np.sum(np.sum(siv_geo_nh, axis=2), axis=1, keepdims=True) 
  
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
selected_volumes=range(len(selected))

for i in range(len(selected)):
    selected_volumes[i]=siv.iloc[0,(selected.iloc[i,0]-min(years))]

PIOMAS=np.loadtxt('/home/josmarti/Data/Observations/PIOMAS_SIV_monthly_1979-2020.txt',usecols=range(1,13))
    
plt.figure(figsize=[12,9])
plt.plot(years, (siv.iloc[0,:]/1e12), label="PPrun")
plt.plot(years[0:41],np.mean((PIOMAS[0:41]),axis=1), label="PIOMAS")
#plt.scatter(selected['Year'], selected_volumes, color='r', s=50)
plt.title('Annual Mean Sea Ice Volume of Control Run \n')
plt.xlabel('Years')
plt.ylabel('SIV')
plt.legend()

"""plt.figure(figsize=[12,9])
plt.plot(years, siv_anoms.mean(axis=0))
plt.scatter(selected['Year'], selected['Mean of Monthly Anomalies'], color='r', s=50)
plt.title('Annual Sea Ice Volume Anomalies of Control Run \n')
plt.xlabel('Years')
plt.ylabel('Anomaly')"""