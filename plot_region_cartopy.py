#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

@author: josmarti
"""
import matplotlib.pyplot as plt
from matplotlib import colorbar, colors
from scipy.io import loadmat
import cartopy.crs as ccrs
import rms_plots as rpl
import numpy as np
import numpy.ma as ma
import nc as nc
import matplotlib as mpl
#from netCDF4 import Dataset

clevslab=np.arange(16)-1 #from rms_plots.py 
r=0
jet3=colors.ListedColormap(loadmat('/home/josmarti/Downloads/cmap_jet3.mat')['cmap'], name='jet3')

####### read data
infile ='/home/josmarti/Data/reg_grid_iceregions.nc'
lsmask=nc.getvar('/home/josmarti/Data/lsmask_cansipsv2_sea.nc', 'LSMASK')
#infile ='/home/josmarti/Data/Observations/had2cis_1x1_198001_202004_sicn.nc'
#infile ='/home/josmarti/Data/sic_monthly_CCCma-CanCM4_NEW_1x1_grid_i20171201.nc'
#infile = '/home/josmarti/Data/Observations/had2cis_1x1_198001_202004_sicn.nc'
#data=Dataset(infile)
#data.set_auto_mask(False)
lon=nc.getvar(infile, 'lon').squeeze()
lat=nc.getvar(infile, 'lat').squeeze()
region=nc.getvar(infile, 'region').squeeze()
#skill=np.load('/home/josmarti/Data/GeoACC_NEW.npy', allow_pickle=True).astype(float)
#mask=nc.getvar('/home/josmarti/Data/1x1_reg_mask.nc', 'region').squeeze()
#mask[mask != ]=0
#region=nc.getvar('/home/josmarti/Data/lsmask_cansipsv2_sea.nc', 'LSMASK') #landsea mask troubleshooting

regionlabs=['0land','1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15OTHER']
region_titles=['Land', 'Central Arctic', 'GIN Seas', 'Barents Sea', 'Kara Sea', 'Laptev Sea', 'East Siberian Sea', 'Chukchi Sea', 'Bering Sea', 'Sea of Okhotsk', 'Beaufort Sea', 'Canadian Archipelago', 'Hudson Bay', 'Baffin Bay', 'Labrador Sea', 'Open Ocean']
if r != 0:
   region=ma.masked_where(region != r, region)

#%%

pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap='acccbar')

####### Plot05
cmap=plt.cm.get_cmap('jet', 16)


fig = plt.figure(figsize=(40, 20))
ax1 = fig.add_subplot(221, projection=ccrs.NorthPolarStereo())
cs1 = ax1.pcolormesh(lon[0:360], lat[120:180], region[120:180,0:360], transform=ccrs.PlateCarree(), cmap=cmap)
#cs1 = ax1.pcolormesh(data['longitude'][0:361], data['latitude'][0:180], np.multiply(data['sic'][0,0,0:180,0:360],np.flip(mask[0:180,0:360])), transform=ccrs.PlateCarree(), cmap='nipy_spectral')
ax1.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.stock_img()

#cb2.set_ticklabels(region_titles)
ax1.set_title("Arctic Regions from Day et al. (2014b)", fontsize=20)
cb=plt.colorbar(cs1)
cb.set_ticks(np.arange(0.4,15.5,0.95))
cb.set_ticklabels(region_titles)
#ax1.set_title('Navy subregions')
#plt.colorbar(cs1,ticks=(range(np.amin(data['region'][:]),np.amax(data['region'][:]))))
#if r == 0:
 #   cbparams=dict(manticks=clevslab, manlabels=region_titles)
  #  rpl.add_cb(ax1,cs1,**pcparams)
