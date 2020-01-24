#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:06:18 2020

@author: Joseph Martin, University of Victoria
"""

from netCDF4 import Dataset
import math as math
import nc as nc
import numpy as np
import array as ar
obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sie')

sim=np.zeros((12,12,32))

#Shape observations
obs=obsin[:,0,0,0]
obs=np.reshape(obs, (12,32))

#Build model array
for months in range(1,13):
	if months<10:
		m="0%s" % months
	else:
		m=str(months)
	for years in range (1979,2011):
		y=str(years)
		var3=nc.getvar(('/home/josmarti/Data/SIE/SIE_monthly_CanCM3_i%s%s.nc' % (y,m)),'sic')
		var3=var3[:,:,0,0] #strips extra dimensions (check which ones)
		var4=nc.getvar(('/home/josmarti/Data/SIE/SIE_monthly_CanCM4_i%s%s.nc' % (y,m)),'sic')
		var4=var4[:,:,0,0]
		var=np.concatenate((var3,var4), axis=1)
		extent=var*2*math.pi*6.371**2 #multiply constant to convert fraction to SIE
		sim[months-1,:,years-1979]=np.mean(extent, axis=1) #Average across ensembles and insert into matrix


#Calculate ACC
ACC=np.zeros((12,12), dtype=np.ndarray)
for init in range(12):
	for lead in range(12):
            ACC[init,lead]=np.corrcoef(sim[init,lead,:],obs[init,:])
print(ACC.shape)

