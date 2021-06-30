#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:21:28 2021

@author: josmarti
"""

import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
import calendar
import rms_plots as rpl
import rms_utils_boot as bt
from scipy.io import loadmat
import Mitch_detrend as md
from scipy.signal import detrend

jet3=colors.ListedColormap(loadmat('/home/josmarti/Downloads/cmap_jet3.mat')['cmap'], name='jet3') #Change filepath (1/2)
mitch=loadmat('/home/josmarti/Downloads/regionalSIE.mat') #Change filepath (2/2)

#regionlabs=['1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15PAN']
region="6ESI" #Choose from regional labels above. 
pstyle="pc" #pc for colour plot, cf for contourf

region_titles=['Central Arctic', 'GIN Seas', 'Barents Sea', 'Kara Sea', 'Laptev Sea', 'East Siberian Sea', 'Chukchi Sea', 'Bering Sea', 'Sea of Okhotsk', 'Beaufort Sea', 'Canadian Archipelago', 'Hudson Bay', 'Baffin Bay', 'Labrador Sea', 'Pan-Arctic']

if region[1].isdigit():
    r=int(region[0:2])
else:
    r=int(region[0])


years=range(1980,2017)
period=max(years) - min(years)
    
obsin=mitch['regionalSIE'][0:456,(r-1)]
obs2mask=obsin[((min(years)-1980)*12):((max(years)+2-1980)*12)]
obs=np.reshape(obs2mask, ((period+2),12)).transpose()

original_obs=obs #used for standard deviation

do_oper_detrend=False #operational detrending 
if do_oper_detrend:
    obs=md.Mitch_detrend(obs)
else:
    for imon in range(12):
        obs[imon,:]=detrend(obs[imon,:])

#######################################################################################################################
auto=np.zeros((12,36))
boot_temp_auto=np.zeros((12,36), dtype=np.ndarray)
boot_auto=np.zeros((12,36), dtype=np.ndarray)

do_Michael=True

if do_Michael:
    for tmon in range(12):
        for lead in range(36):
            imon=tmon-lead%12-1
            yearoffset=lead/12    
            if imon<0:
                yearoffset=yearoffset+1
                imon=imon+12
            #print tmon,lead,imon,yearoffset
            auto[tmon,lead]=np.corrcoef(obs[imon,(3-yearoffset):(38-yearoffset)],obs[tmon,3::])[1,0]
            boot_temp_auto[tmon,lead]=bt.calc_corr_boot(obs[imon,(3-yearoffset):(38-yearoffset)],obs[tmon,3::], 1000)
            boot_auto[tmon,lead]=bt.calc_boot_stats(boot_temp_auto[tmon,lead],sides=1,pval_threshold=0.05)[1]
            
"""
for i in range(12):
    for j in range(36):
        if (i-j-1)<0:
            auto[i,j]=np.corrcoef(obs[(i-((j % 12)+1)),0:-((j/12)+1)],obs[i,((j/12)+1)::])[1,0]
            boot_temp_auto[i,j]=bt.calc_corr_boot(obs[(i-((j % 12)+1)),0:-((j/12)+1)],obs[i,((j/12)+1)::], 1000)
            boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]
        else:
            if (j/12)==0:
                auto[i,j]=np.corrcoef(obs[(i-((j % 12)+1)),:],obs[i,:])[1,0]
                boot_temp_auto[i,j]=bt.calc_corr_boot(obs[(i-((j % 12)+1)),:],obs[i,:], 1000)
                boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]
            else:
                auto[i,j]=np.corrcoef(obs[(i-((j % 12)+1)),0:-(j/12)],obs[i,(j/12)::])[1,0]
                boot_temp_auto[i,j]=bt.calc_corr_boot(obs[(i-((j % 12)+1)),0:-(j/12)],obs[i,(j/12)::], 1000)
                boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]"""

#######################################################################################################################
                
if r!=0: #Masks out for STD
    if (len(np.std(original_obs, axis=1)[(np.std(original_obs, axis=1))< 0.03e12])!=0):
        print ("STD mask in effect")
        for i in range(len(auto)):
            if np.std(original_obs, axis=1)[i] < 0.03e12:
                boot_auto[i,:]=-1  # 0 will be seen as significant
                auto[i,:]=0
#%%
                
fig, ax = plt.subplots(figsize=(4,9))
if pstyle == "pc":
    #pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar') #CCCMA
    pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)  #Mitch's
    pc=rpl.add_pc(ax,range(13),range(37),auto.transpose(),**pcparams)
elif pstyle == "cf":
    pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3,latlon=False)
    pc=rpl.add_cf(ax,range(1,13),range(1,37),auto.transpose(),**pcparams)
for month in range(len(boot_auto[:,0])):
            for lead in range(len(boot_auto[0,:])):
                if boot_auto[month,lead]>0:   #MUST BE CHANGED BACK to >=!!!!!!
                    plt.scatter((month+0.5),(lead+0.5), color = 'black', s=20)
                    
plt.title("%s Obs AutoCorr" % region_titles[(r-1)], fontsize=16)
plt.colorbar(pc)
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
ax.set_yticks(np.arange(36)+0.5)    
ax.set_yticklabels(range(36))  
plt.xlabel("Predicted Month")    
plt.ylabel("Lead (Months)")      
plt.show()   