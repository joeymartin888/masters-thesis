#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 21:14:57 2021

@author: josmarti
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import cartopy.crs as ccrs
import rms_plots as rpl
import nc as nc
import numpy as np
import re
import calendar
import rms_utils_boot as bt
import CCMA_plot
from scipy.io import loadmat
import math as math

years=range(1979,2011)

monthly_matrix=np.zeros(((10*(max(years)+1-min(years))),12))
model=np.zeros(((10*(max(years)+1-min(years))),12))
control=np.zeros((12,max(years)+1-min(years)))
ACC=np.zeros(12, dtype=np.ndarray)
for y in years:
    for i in range(1,11):
        monthly_matrix[((y-min(years))*10+(i-1)),:]=nc.getvar(('/home/josmarti/Data/APPOSITE/pred/Jul/SIE/SIE_Omon_pred_%iJul%i.nc' % (y,i)), 'sic').squeeze()*2*math.pi*6.371**2

for k in years:
    control[:,k-min(years)]=nc.getvar(('/home/josmarti/Data/APPOSITE/ctrl/sie_Omon_ctrl_%i.nc' % k), 'sic').squeeze()*2*math.pi*6.371**2

clim_mean=np.mean(control, axis=1)
#%%
for y in range(len(years)):
    for e in range(10):
        model[(y*10+e),:]=np.mean(np.delete(monthly_matrix[y*10:y*10+10,:], e, axis=0), axis=0)

#%%    
for j in range(12):
    ACC[j]=np.corrcoef(model[:,j],monthly_matrix[:,j])[1,0]

plt.plot(ACC)