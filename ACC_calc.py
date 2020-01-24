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
print(obs.shape)
print(np.mean(obs))
#print(np.mean(var*2*math.pi*6.371**2))
#print(var.shape)
#print(obs)

model="CanCM4"
init="198812"
#Model into ACC matrix format
var=nc.getvar(('/home/josmarti/Data/SIE/SIE_monthly_%s_i%s.nc' % (model,init)),'sic')
#print(var)

#Observations into years
for i in range(32):
    index_pos=range(i*12,i*12+12)
    obs=obsin[index_pos]
    #for j in range(12):
        #print(j)
        
#varem=np.mean(var, axis=1) #Average across ensembles
#print(varem)
#ACC=np.corrcoef(sim,obs)

print(ACC)

