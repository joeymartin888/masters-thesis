#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

@author: josmarti
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import rms_plots as rpl
import numpy as np
import numpy.ma as ma
import nc as nc
from netCDF4 import Dataset

clevslab=np.arange(18)-1 #from rms_plots.py 
r=0

####### read data
infile = '/home/josmarti/Data/1x1_reg_mask.nc'
data=Dataset(infile)
data.set_auto_mask(False)

regionlabs=['0land','1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15OTHER']
if r != 0:
    region[region != r]=np.nan


####### Plot05
fig = plt.figure(figsize=(40, 20))
ax1 = fig.add_subplot(221, projection=ccrs.NorthPolarStereo())
cs1 = ax1.pcolormesh(data['lon'][0:361], data['lat'][120:180], data['region'][0,0,120:180,0:361], transform=ccrs.PlateCarree(), cmap='nipy_spectral')
ax1.set_extent([-180, 180, 48, 90], crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.stock_img()
#ax1.set_title('Navy subregions')
#plt.colorbar(cs1,ticks=(range(np.amin(data['region'][:]),np.amax(data['region'][:]))))
if r == 0:
    cbparams=dict(manticks=clevslab, manlabels=regionlabs)
    rpl.add_cb(ax1,cs1,**cbparams)

