#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 12:51:48 2021

@author: josmarti
"""

import rms_plots as rpl
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
import calendar


version=raw_input("Version?\n")
skill=np.load('/home/josmarti/Data/GeoACC_%s.npy' % version, allow_pickle=True).astype(float)
ACC=np.nanmean(np.nanmean(skill, axis=2),axis=2)
rpl.register_rms_cmaps()
fig, ax = plt.subplots()
pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
pc=rpl.add_pc(ax,range(13),range(13),ACC,**pcparams)
plt.colorbar(pc)
plt.title("ACC Means for %s" % version)
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
ax.set_yticklabels(['0','1','2','3','4','5','6','7','8','9','10','11'])  
plt.xlabel("Predicted Month")    
plt.ylabel("Lead (Months)")      
plt.show()