#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:06:18 2020

@author: Joseph Martin, University of Victoria 
with contributions from Dr. Michael Sigmond, Canadian Centre for Climate Modelling and Analysis

"""

from netCDF4 import Dataset
import math as math
import nc as nc
import numpy as np
import array as ar
from matplotlib import colors
import matplotlib.pyplot as plt
import calendar
import rms_plots as rpl
import CCMA_plot
from scipy import signal
import rms_utils_boot as bt
from scipy.io import loadmat
import Mitch_detrend as md

#regionlabs=['0land','1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15OTHER']
region="0NONE" 

if region[1].isdigit():
    r=int(region[0:2])
else:
    r=int(region[0])

#Select metric
metric="SIE"

#Select observation data set
Data="Had2CIS"

#obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', str.lower(metric)).squeeze()

#Select Year Range
years=range(1980,2019)
period=max(years) - min(years)

sim=np.zeros((12,12,len(years)))

#Shape observations
if r !=0 :
    """obsingeo=nc.getvar('/home/josmarti/Data/Observations/had2cis_1x1_198001_202004_sicn.nc', 'SICN').squeeze()
    #obsingeo[obsingeo.mask==True]=1 #CHECK WITH MICHAEL!!!!!!
    if metric=="SIE":
        obsingeo[obsingeo <= 0.15] = 0
        obsingeo[obsingeo > 0.15] = 1
    mask=nc.getvar('/home/josmarti/Data/1x1_reg_mask.nc', 'region').squeeze()
    gridpoint=nc.getvar('/home/josmarti/Data/gridpoint_1x1.nc','areacella')
    mask[mask != r]=0
    mask[mask == r]=1
    obs1=np.multiply(obsingeo.astype(float),mask.astype(float)) #changes to dtype=float64 to avoid run time error
    #obs2=obs1 #Select Northern Hemisphere
    obsnh=np.multiply(obs1,gridpoint.astype(float))
    obstemp=np.sum(np.sum((obsnh/(1e6*2.8**2)), axis=1), axis=1) #CHANGED FROM MEAN"""
    obsin=nc.getvar(('/home/josmarti/Data/Observations/Observed_SIE_%i.nc' % r), 'SICN').squeeze()*2*math.pi*6.371**2
    obs2mask=obsin[((min(years)-1980)*12):((max(years)+2-1980)*12)]
    obs=np.reshape(obs2mask, ((period+2),12)).transpose()
else:
    obsingeo=np.delete(nc.getvar('/home/josmarti/Data/Observations/had2cis_128_64_195901_202004_sic.nc', 'SICN').squeeze(), 128, 2)
    gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')
    obs1=np.multiply(obsingeo,gridpoint)
    obsnh=np.delete(obs1, range(32), axis=1)
    obstemp=np.mean(np.mean(obsnh, axis=1), axis=1)
    obs2mask=obstemp[((min(years)-1959)*12):((max(years)+2-1959)*12)] 
    if Data == "Had2CIS":
        if metric == "SIA":
            obs=np.reshape(obs2mask, ((period+2),12)).transpose()
        elif metric == "SIE":
            obsin=nc.getvar('/home/josmarti/Data/Observations/had2cisSIE.nc', 'SICN').squeeze()[((min(years)-1959)*12):((max(years)+2-1959)*12)] 
            obs=np.reshape(obsin, ((period+2),12)).transpose()
    elif Data == "NSIDC":
        obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', str.lower(metric)).squeeze()
        obs=np.delete((np.reshape(obsin, ((period+3),12)).transpose()),range(min(years)-1979),1)
        


#%%

#Option for linear detrending
detrend=True

#Display persistence forecast
show_persistence=False

#Mask Low STD
std_mask=True

if std_mask:
    mask=nc.getvar('/home/josmarti/Data/1x1_reg_mask.nc', 'region').squeeze()
    gridpoint=nc.getvar('/home/josmarti/Data/gridpoint_1x1.nc','areacella')
    mask[mask != r]=0
    mask[mask == r]=1
    area=np.sum(np.multiply(gridpoint,mask))

#Set figure titles for thesis
thesis_figures=True
titles=[]

#Mitch's colors
jet3=colors.ListedColormap(loadmat('/home/josmarti/Downloads/cmap_jet3.mat')['cmap'], name='jet3') 

#Used to seperate model outputs
single=False


#Select CANSIPS v1 or v2
CANSIPS="v1"
#CANSIPSes=["v1", "v2"]
#CANSIPS="v2"

#Select OLD or NEW
#version="NEW"
#version="OLD"
versions=["OLD", "NEW"]
#versions=["OLD/1x1_trial", "NEW/1x1_trial"]

#Set plot style
pstyle = "pc"
if pstyle != "cf" and pstyle != "pc":
    print("Select a plot style (cf or pc)")
    
#Build persistance model
persistence=np.zeros(sim.shape)
obs_pers=np.zeros(obs.shape)
obs_pers[0:11,:]=obs[0:11,:]
#if r!=0:
if Data == "Had2CIS":
    obs_pers[11,0]=obsin[(min(years)-1-1959)*12+11]
elif Data == "NSIDC":
    obs_pers[11,0]=obsin[(min(years)-1980)*12+11]
#else:    
 #   obs_pers[11,0]=obsin[(min(years)-1980)*12+11]
obs_pers[11,1::]=obs[11,0:-1]
obs_mean=np.mean(obs_pers, axis=1, keepdims=True)
obs_anom=obs_pers-obs_mean
for init in range(12):
	for target in range(12):
            persistence[init,target,:]=obs_mean[target]+obs_anom[init-1,0:-1] 

#%%
            
original_obs=obs #used for standard deviation

if detrend: #sim detrended below
    persistence=signal.detrend(persistence)
    obs=signal.detrend(obs)
    #obs=md.Mitch_detrend(obs)

diff=0  
#%%    
for version in versions:
#for CANSIPS in CANSIPSes:
#for smark in range(2):    
    #Select SIE or SIA

        
    """if version=="NEW":
        CANSIPS="v2"
    #for overall comparison"""
    
    
    #Build model array
    for months in range(1,13):
    	if months<10:
    		m="0%s" % months
    	else:
    		m=str(months)
    	for year in years:
                y=str(year)
                if CANSIPS=="v1":
                    if r !=0 :
                        var3=nc.getvar(('/home/josmarti/Data/%s/%s/%s/%s_monthly_%s_CanCM3_i%s%s.nc' % (version,metric,region,metric,region,y,m)),'sic').squeeze()
                    else:
                        var3=nc.getvar(('/home/josmarti/Data/%s/%s/%s_monthly_CanCM3_i%s%s.nc' % (version,metric,metric,y,m)),'sic').squeeze()#.transpose()
                if CANSIPS=="v2":
                    if r!=0:
                        var3=nc.getvar(('/home/josmarti/Data/GEM-NEMO/%s/%s/%s_monthly_%s_GEM_NEMO_i%s%s.nc' % (metric,region,metric,region,y,m)),'sic').squeeze()
                    else:
                        var3=nc.getvar(('/home/josmarti/Data/GEM-NEMO/%s/%s_monthly_GEM-NEMO_i%s%s.nc' % (metric,metric,y,m)),'sic').squeeze()
                #if r!=0:
                #    var4=nc.getvar(('/home/josmarti/Data/%s/%s/%s/%s_monthly_%s_CanCM4_i%s%s.nc' % (version,metric,region,metric,region,y,m)),'sic').squeeze()
                #else:
                if r !=0 :
                    var4=nc.getvar(('/home/josmarti/Data/%s/%s/%s/%s_monthly_%s_CanCM4_i%s%s.nc' % (version,metric,region,metric,region,y,m)),'sic').squeeze()
                else:
                    var4=nc.getvar(('/home/josmarti/Data/%s/%s/%s_monthly_CanCM4_i%s%s.nc' % (version,metric,metric,y,m)),'sic').squeeze()#.transpose()
                var=np.concatenate((var3,var4), axis=1)
                if single:
                    print("Single Model")
                    if smark==0:
                        var=var3
                    elif smark==1:
                        var=var4
                if metric=="SIA":
                    sim[months-1,:,year-min(years)]=np.mean(var, axis=1) #Average across ensembles and insert into matrix
                if metric=="SIE":
                    extent=var*2*math.pi*6.371**2 #multiply constant to convert fraction to SIE
                    #avex=np.mean(extent, axis=1)
                    sim[months-1,:,year-min(years)]=np.mean(extent, axis=1) #Average across ensembles and insert into matrix
#%%
    
    for i in range(12): #line up all target months
        sim[i,:,:]=np.roll(sim[i,:,:],i,axis=0)
    
   

                    
        
    if detrend:
        #sim=signal.detrend(sim)         #Linear detrending, obs detrended above 
        sim=md.Mitch_detrend(sim)
                
        """for i in range(12):                #Polynomial detrending
            for j in range(12):
                t_sim=np.arange(len(sim[i,j,:]))
                d_sim=np.polyfit(t_sim, sim[i,j,:], 1)
                sim[i,j,:]=sim[i,j,:]-d_sim[0]*t_sim
            t_obs=np.arange(len(obs[i,:]))
            d_obs=np.polyfit(t_obs, obs[i,:], 1)
            obs[i,:]=obs[i,:]-d_obs[0]*t_obs   """ 
       
    #print(persistence)
    
    
    #Calculate ACC
    ACC=np.zeros((12,12), dtype=np.ndarray)
    std_dev=np.zeros((12,12), dtype=np.ndarray)
    boot_temp=np.zeros((12,12), dtype=np.ndarray)
    boot=np.zeros((12,12), dtype=np.ndarray)
    ACC_pers=np.zeros((12,12), dtype=np.ndarray)
    boot_temp_pers=np.zeros((12,12), dtype=np.ndarray)
    boot_pers=np.zeros((12,12), dtype=np.ndarray)
    for init in range(12):
    	for target in range(12):
                if init<=target: 
                    ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,0:-1])[1,0]
                    std_dev[init,target]=np.std(original_obs[target,0:-1])
                    boot_temp[init,target]=bt.calc_corr_boot(sim[init,target,:],obs[target,0:-1],1000)
                    boot[init,target]=bt.calc_boot_stats(boot_temp[init,target],sides=1,pval_threshold=0.05)[1]
                    ACC_pers[init,target]=np.corrcoef(persistence[init,target,:],obs[target,0:-1])[1,0]
                    boot_temp_pers[init,target]=bt.calc_corr_boot(persistence[init,target,:],obs[target,0:-1],1000)
                    boot_pers[init,target]=bt.calc_boot_stats(boot_temp_pers[init,target],sides=2,pval_threshold=0.05)[1]
                else: #rolls observations to realign years
                    ACC[init,target]=np.corrcoef(sim[init,target,:],obs[target,1::])[1,0]
                    std_dev[init,target]=np.std(original_obs[target,1::])
                    boot_temp[init,target]=bt.calc_corr_boot(sim[init,target,:],obs[target,1::],1000)
                    boot[init,target]=bt.calc_boot_stats(boot_temp[init,target],sides=1,pval_threshold=0.05)[1]
                    ACC_pers[init,target]=np.corrcoef(persistence[init,target,:],obs[target,1::])[1,0]
                    boot_temp_pers[init,target]=bt.calc_corr_boot(persistence[init,target,:],obs[target,1::],1000)
                    boot_pers[init,target]=bt.calc_boot_stats(boot_temp_pers[init,target],sides=2,pval_threshold=0.05)[1]
    
    #print(ACC_pers)
    
    
    
          

              
                
    #boot[init,target]=bt.calc_corr_boot(sim,obs,1000)
    
   
    
    
    #print(len(boot[5,3]))
    """nt=np.shape(sim)[0]  
    nboot=np.array(1000,dtype=int)
    dims=np.concatenate(([nboot],np.array(np.shape(sim)[1::])))
    print(dims)"""

    ACC = np.vstack(ACC[:, :]).astype(np.float) #ACC is an object and pcolor needs floats
    ACC_pers = np.vstack(ACC_pers[:, :]).astype(np.float) #ACC_pers is an object and pcolor needs floats
    
    
    ACC2=np.zeros((12,12))
    boot_temp2=np.zeros((12,12), dtype=np.ndarray)
    boot2=np.zeros((12,12))
    ACC2_pers=np.zeros((12,12))
    boot_temp2_pers=np.zeros((12,12), dtype=np.ndarray)
    boot2_pers=np.zeros((12,12))

    for init in range(12):
        for target in range(12):
            if target>=init: #same year
                ilag=target-init
            else: #next year 
                ilag=target-init+12
            ACC2[ilag,target]=ACC[init,target]
            boot_temp2[ilag,target]=boot_temp[init,target]
            boot2[ilag,target]=boot[init,target]
            ACC2_pers[ilag,target]=ACC_pers[init,target]
            boot_temp2_pers[ilag,target]=boot_temp_pers[init,target]
            boot2_pers[ilag,target]=boot_pers[init,target]
            
            
    if diff==0:
        old_ACC=ACC2
        old_boot=boot_temp2
        diff+=1
    
    if diff==1:
        new_ACC=ACC2
        new_boot=boot_temp2
    
    if r!=0:
        if (len(np.std(original_obs, axis=1)[(np.std(original_obs, axis=1)*1e15)/np.sum(area) < 0.8])!=0 and std_mask):
            print ("STD mask in effect")
            for i in range(len(ACC2)):
                if ((np.std(original_obs, axis=1)[i]*1e15)/np.sum(area)) < 0.8:
                    boot2_pers[:,i]=-1  # 0 will be seen as significant
                    boot2[:,i]=-1 # 0 will be seen as significant
                    ACC2_pers[:,i]=0
                    ACC2[:,i]=0
    
    if show_persistence:    
        fig, ax = plt.subplots()
        
        pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)
        #pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
        pc=rpl.add_pc(ax,range(13),range(13),ACC2_pers,**pcparams)
        plt.colorbar(pc)
        for init in range(len(boot2_pers[:,0])):
            for target in range(len(boot2_pers[0,:])):
                if boot2_pers[init,target]>=0:
                    plt.scatter((target+0.5),(init+0.5), color = 'black', s=20, marker='^')
        if detrend:
            plt.title("ACC for Detrended Persistence from %i to %i" % (min(years),(max(years)+1)))
        else:
            plt.title("ACC for Persistence from %i to %i" % (min(years),(max(years)+1)))
        ax.invert_xaxis
        ax.invert_yaxis
        ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
        ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
        ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
        ax.set_yticklabels(['0','1','2','3','4','5','6','7','8','9','10','11'])  
        plt.xlabel("Predicted Month")    
        plt.ylabel("Lead (Months)")      
        plt.show()   

#%%
    print("%s" % Data)
    fig, ax = plt.subplots()
    if pstyle == "pc":
        #pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)
        pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
        pc=rpl.add_pc(ax,range(13),range(13),ACC2,**pcparams)
        #ss=rpl.add_sc(ax,range(13),range(13),boot2)
    elif pstyle == "cf":
        pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar', latlon=False)
        pc=rpl.add_cf(ax,range(12),range(12),ACC2,**pcparams)
    plt.colorbar(pc)
    for init in range(len(boot2[:,0])):
        for target in range(len(boot2[0,:])):
            if boot2[init,target]==np.nan:
                plt.scatter((target+0.5),(init+0.5), color = 'red', s=50, marker='x')
                print("done")
            if boot2[init,target]>0:
                if ACC2_pers[init,target]>=ACC2[init,target]:
                    plt.scatter((target+0.5),(init+0.5), color = 'black', s=20)
                else:
                    plt.scatter((target+0.5),(init+0.5), color = 'black', s=40, marker='^')
            """else:
                if ACC2[init,target]>ACC2_pers[init,target]:
                    plt.scatter((target+0.5),(init+0.5), color = 'blue', s=20, marker='s')"""
    if detrend:
        if r!=0:
            plt.title("%s %s %s ACC of %s (detrended) from %i to %i" % (CANSIPS, str.capitalize(version), metric, region, min(years),(max(years))))
        else:
            plt.title("%s %s %s ACC of detrended forecasts from %i to %i" % (CANSIPS, str.capitalize(version), metric, min(years),(max(years))))
    else:
        if r!=0:
            plt.title("%s %s %s ACC of %s from %i to %i" % (CANSIPS, str.capitalize(version), metric, region, min(years),(max(years))))
        else:
            plt.title("%s %s %s ACC of forecasts from %i to %i" % (CANSIPS, str.capitalize(version), metric, min(years),(max(years))))
    if single:
        if smark==0:
            plt.title("GEM-NEMO %s ACC of forecasts from %i to %i" % (metric, min(years),(max(years))))
        elif smark==1:
            if r!=0:
                plt.title("%s Operational Skill" % region, fontsize=20)
            else:
                plt.title("CanCM4 Operational Skill", fontsize=20)
    if thesis_figures:
        print("%s" % metric)
        if CANSIPS=="v1":
            if version=="OLD":
                plt.title("CanSIPSv1", fontsize=20)
                titles.append("CanSIPSv1")
            elif version=="NEW":
                plt.title("CanSIPSv1b", fontsize=20)
                titles.append("CanSIPSv1b")
        elif CANSIPS=="v2":
                plt.title("CanSIPSv2", fontsize=20)
                titles.append("CanSIPSv2")
    ax.invert_xaxis
    ax.invert_yaxis
    ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
    ax.set_xticklabels(calendar.month_name[1:13], rotation=90, fontsize=12)    
    ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
    ax.set_yticklabels(['0','1','2','3','4','5','6','7','8','9','10','11'], fontsize=12)    
    plt.ylabel("Lead (months)", fontsize=15)
    plt.xlabel("Target Month", fontsize=15)      
    plt.show()
    
    if single:
        if smark == 1:
            np.save('temp/%s' % region, ACC2)
    
#%%     
   
    
    """
    k=pd.DataFrame(ACC2,index=calendar.month_name[1:13], columns=range(1,13))
    fig, ax = plt.subplots()
    d=ax.pcolor(k)
    plt.colorbar(d)
    plt.title("%s %s ACC of forecasts from %i to %i" % (str.capitalize(version), metric, min(years),(max(years)+1)))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    #ax.xaxis.set_minor_locator(ticker.FixedLocator((range(1,13)+0.5)))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(str(range(1,13))))
    ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)
    ax.set_yticklabels(k.index)
    plt.xlabel("Lead Time")
    plt.ylabel("Initialization month")
    plt.show() """

#%%  
ACC=new_ACC-old_ACC
boot_temp=new_boot-old_boot
for init in range(12):
        for target in range(12):
            boot[init,target]=bt.calc_boot_stats(boot_temp[init,target],sides=1,pval_threshold=0.05)[1]

if r!=0:
        if (len(np.std(obs, axis=1)[((np.std(obs, axis=1)*1e15)/np.sum(area)) < 0.8])!=0 and std_mask):
            for i in range(len(ACC)):
                if ((np.std(obs, axis=1)[i]*1e15)/np.sum(area)) < 0.8:
                    boot[:,i]=-1 # 0 will be seen as significant


fig, ax = plt.subplots()
if pstyle == "pc":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar')
    pc=rpl.add_pc(ax,range(13),range(13),ACC,**pcparams)
elif pstyle == "cf":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar',latlon=False)
    pc=rpl.add_cf(ax,range(1,13),range(1,13),ACC,**pcparams)
plt.colorbar(pc)
for init in range(len(boot[:,0])):
    for target in range(len(boot[0,:])):
        if boot[init,target]>0:
            plt.scatter((target+0.5),(init+0.5), color = 'black', s=40, marker='^')
if detrend:
    plt.title("Difference in ACC of detrended forecasts from %i to %i" % (min(years),(max(years))))
else:
    plt.title("Difference in ACC of forecasts from %i to %i" % (min(years),(max(years))))
if thesis_figures:
    plt.title("%s - %s" % (titles[1],titles[0]), fontsize=20)
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90, fontsize=12)    
ax.set_yticks(np.arange(len(calendar.month_name[1:13]))+0.5)    
ax.set_yticklabels(['0','1','2','3','4','5','6','7','8','9','10','11'], fontsize=12)  
plt.ylabel("Lead (months)", fontsize=15)
plt.xlabel("Target Month", fontsize=15)    
plt.show()
