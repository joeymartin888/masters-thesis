#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 13:56:25 2020

@author: josmarti
"""

import sys
sys.path.append('/HOME/rms/PMODS')
import os as os
import nc as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

os.system('cdo selmon,3 ~/Data/Observations/had2cis_1x1_198001_202004_sicn.nc ./tmp1.nc')
os.system('cdo remapbil,1x1_reg_mask.nc tmp1.nc tmp2.nc') #regrid
os.system('cdo -gtc,0.15 tmp2.nc tmp3.nc') # SIE mask: 1's where sic>0.15, 0 elsewhere

os.system('cdo -setctomiss,0 -eqc,1 1x1_reg_mask.nc 1centralarctic.nc') # create central arctic mask

os.system('cdo mul 1centralarctic.nc tmp3.nc tmp4.nc') # mask out everything outside central Arctic

os.system('cdo -fldmean -sellonlatbox,0,360,0,90 tmp4.nc tmp5.nc') # Average fraction per km^2 in NH 

SIE_im03_ca=nc.getvar('tmp5.nc','SICN').squeeze()*2*np.pi*6371.*6371./1e6

plt.plot(signal.detrend(SIE_im03_ca))
print np.std(signal.detrend(SIE_im03_ca))