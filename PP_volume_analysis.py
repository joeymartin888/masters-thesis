#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

@author: josmarti
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import rms_plots as rpl
import nc as nc
import numpy as np
import re
#import cdo; c=cdo.Cdo()

syear=1980; eyear=1995

obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sia').squeeze()
start=(1990-1979)*12; end=start+12
mean_obs=np.mean(np.reshape(obsin[((1980-1979)*12):((1995-1979)*12)],(15,12)).transpose(), axis=1)
gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')

PIOMAS=np.loadtxt('/home/josmarti/Data/Observations/PIOMAS_SIV_monthly_1979-2020.txt',usecols=range(1,13))
mean_obs=np.mean(PIOMAS[(syear-1979):(eyear-1979)], axis=0)

years=[2312, 2435, 2413, 2356, 2387, 2332]
ensembles=["%03d" % i for i in range(1,13)]
#ensembles=["001"]
run_length=1 #years
daily_matrix=np.zeros((6,12,1095))
days=365*3

for e in ensembles:
    for year in years:
        var1=np.multiply(nc.getvar(('/home/josmarti/Data/PPrun/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year),str(year))), 'sic'), gridpoint/913)
        var2=np.multiply(nc.getvar(('/home/josmarti/Data/PPrun/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+1),str(year+1))), 'sic'), gridpoint/913)
        var3=np.multiply(nc.getvar(('/home/josmarti/Data/PPrun/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+2),str(year+2))), 'sic'), gridpoint/913)
        var=nc.getvar('/home/josmarti/Data/PPrun_Jan/sc_dhfp1e_e%s_i%s_m01_daily_sicn.nc4' % (e,str(year)),'sicn')
        var1nh=np.delete(var1, range(32), axis=1)
        var2nh=np.delete(var2, range(32), axis=1)
        var3nh=np.delete(var3, range(32), axis=1)
        geosum1=np.sum(np.sum(var1nh, axis=1), axis=1)
        geosum2=np.sum(np.sum(var2nh, axis=1), axis=1)
        geosum3=np.sum(np.sum(var3nh, axis=1), axis=1)
        daily_matrix[years.index(year),ensembles.index(e),:]=np.sum(np.sum(var[:,32:63,:], axis=1), axis=1)
        if year==2413 and e=='011': #2415 e11 seems to have an extra timestep
            geosum3=np.delete(geosum3, 12)
        #run=np.concatenate([geosum1, geosum2, geosum3])
        run=geosum1
        if len(run)>36:
            print("%s, %s" % (str(year), e))
        months=range(1,13)
        plt.plot(months, geosum1, label=('%se%s' % (str(year),e)))
        if run_length==2:
            plt.plot(months, geosum2, label=('%se%s' % (str(year+1),e)))
        if run_length==3:
            plt.plot(months, geosum2, label=('%se%s' % (str(year+1),e)))
            plt.plot(months, geosum3, label=('%se%s' % (str(year+2),e)))
#plt.plot(range(1,13), (np.mean(run/obsin[start:end])*obsin[start:end]), label='observations')
plt.plot(months, (1e12*mean_obs), label='mean_obs',  lw=2.5, color='black')
plt.title("%i year run of %i ensembles (observations: %i-%i)" % (run_length,len(ensembles),syear,eyear))
if (len(years)*run_length*len(ensembles))<12:
    plt.legend()
#plt.legend()

#%%
"""    
fig=plt.figure(figsize=(10,10))

for y in years:
    fig.add_subplot(3,2,(years.index(y)+1), title="Std of %s" % str(y))
    for i in range(12):
        #plt.plot(range(days), daily_matrix[years.index(y),i,0:days], label='%s' % ensembles[i])
        plt.plot(range(days), np.std(daily_matrix[years.index(y),:,0:days], axis=0))
    #plt.title("First %i daily SIA values of %i for all 12 ensembles" % (days,y))
    plt.legend()
fig.tight_layout()"""