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
import nc as nc

clevslab=np.arange(17)-1 #from rms_plots.py 

####### read data
infile = '/home/josmarti/Data/iceregions_128_64.nc'
lon=nc.getvar(infile,'lon')
lat=nc.getvar(infile,'lat')
region=nc.getvar(infile,'REG').squeeze()
regionlabs=['-1ocean','0land','1ARC','2BER','3STL','4BAF','5GRE','6BAR','7KAR','8LAP','9ESI','10CHU','11BEA','12CAN','13HUD','14OKH']

####### Plot
fig = plt.figure(figsize=(40, 20))

ax1 = fig.add_subplot(221, projection=ccrs.NorthPolarStereo())
cs1 = ax1.contourf(lon, lat, region, 20, transform=ccrs.PlateCarree())
ax1.set_extent([-180, 180, 48, 90], crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.stock_img()
ax1.set_title('Navy subregions')

cbparams=dict(manticks=clevslab, manlabels=regionlabs)
rpl.add_cb(ax1,cs1,**cbparams)




