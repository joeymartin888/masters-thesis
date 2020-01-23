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
obs=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', 'sie')
var=nc.getvar('SIE_monthly_CanCM3_i197902.nc','sic')
print(obs)
print(np.mean(obs))
print(np.mean(var*2*math.pi*6.371**2))
#print(var)
#print(obs)

f=4
c=2
a=6

"""
f=ar.array('f',[1,2,3])
c=ar.array('f',[2,3,4])
a=ar.array('f',[3,4,5])
"""

ACCtop=np.mean(np.multiply(np.subtract(f,c),np.subtract(a,c)))
ACCbottom=np.sqrt(np.mean(np.multiply(np.square(np.subtract(f,c)),np.square(np.subtract(a,c)))))
ACC=np.divide(ACCtop,ACCbottom)

print(ACC)
