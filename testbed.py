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

years=range(1980,2010)
r=13

obsingeo=np.delete(nc.getvar('/home/josmarti/Data/Observations/had2cis_128_64_195901_202004_sic.nc', 'SICN').squeeze(), 128, 2)
obs2mask=obsingeo[((min(years)-1959)*12):((max(years)+2-1959)*12),:,:]
mask=nc.getvar('/home/josmarti/Data/iceregions_128_64.nc', 'REG').squeeze()
gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')
for m in range(len(mask)):
    for n in range(len(mask[0])):
        if mask[m,n] != r:
            mask[m,n]=0
        elif mask[m,n] == r:
            mask[m,n]=1
obs1=np.multiply(np.multiply(mask,obs2mask),gridpoint)
obsnh=np.delete(obs1, range(32), axis=1)
obstemp=np.mean(np.mean(obsnh, axis=1), axis=1)
obs=np.reshape(obstemp, (31,12)).transpose()
#var=nc.getvar('/home/josmarti/Projects/masters-thesis/tmp3.nc', 'sic')







