#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

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
#import cdo; c=cdo.Cdo()

auto_corr=True
show_timeseries=True
thesis_figures=True

region_titles=['Pan-Arctic', 'Central Arctic', 'GIN Seas', 'Barents Sea', 'Kara Sea', 'Laptev Sea', 'East Siberian Sea', 'Chukchi Sea', 'Bering Sea', 'Sea of Okhotsk', 'Beaufort Sea', 'Canadian Archipelago', 'Hudson Bay', 'Baffin Bay', 'Labrador Sea']
#regionlabs=['0land','1ARC','2GIN','3BAR','4KAR','5LAP','6ESI','7CHU','8BER','9OKH','10BEA','11CAN','12HUD','13BAF','14LAB','15OTHER']
region="0NONE"
if region[1].isdigit():
    r=int(region[0:2])
else:
    r=int(region[0])
    
std_mask=True

if std_mask:
    mask=nc.getvar('/home/josmarti/Data/1x1_reg_mask.nc', 'region').squeeze()
    gridpoint=nc.getvar('/home/josmarti/Data/gridpoint_1x1.nc','areacella')
    mask[mask != r]=0
    mask[mask == r]=1
    area=np.sum(np.multiply(gridpoint,mask))
#%%
cyears=range(2300,2449)
control=np.zeros((12,len(cyears)))

if r!=0:
    for year in cyears:
        control[:,cyears.index(year)]=nc.getvar(('/home/josmarti/Data/CanCM4_control/SIE/%s/sc_dhfp1e_e001_%s_%i_SIE.nc' % (region,region,year)), 'sicn').squeeze()*2*math.pi*6.371**2
else:    
    for year in cyears:
        control[:,cyears.index(year)]=nc.getvar(('/home/josmarti/Data/CanCM4_control/SIE/sc_dhfp1e_e001_%i_SIE.nc' % year), 'sicn').squeeze()*2*math.pi*6.371**2

#%%

data=loadmat('/home/josmarti/Downloads/SIE_16.mat')

pstyle="pc"

metric='SIV'

syear=1979; eyear=2016; 

obsin=np.zeros((12,40)) 

if metric=='SIE':
    for m in ["%02d" % i for i in range(1,13)]:
            obsin[int(m)-1,:]=np.loadtxt('/home/josmarti/Data/Observations/NSIDC/north/N_%s_extent_v3.0.csv' % m,usecols=4, skiprows=1, delimiter=',', max_rows=40)
    obsin[obsin == -9999] = np.nan
if metric=='SIA':
    for m in ["%02d" % i for i in range(1,13)]:
            obsin[int(m)-1,:]=np.loadtxt('/home/josmarti/Data/Observations/NSIDC/north/N_%s_extent_v3.0.csv' % m,usecols=5, skiprows=1, delimiter=',', max_rows=40)
    obsin[obsin == -9999] = np.nan
mean_obs=np.nanmean(obsin[:,(syear-1979):(eyear-1979)], axis=1)*1e12
annual_mean_obs=np.nanmean(obsin[:,(syear-1979):(eyear-1979)], axis=0)*1e12
#gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')

PIOMAS=np.loadtxt('/home/josmarti/Data/Observations/PIOMAS_SIV_monthly_1979-2020.txt',usecols=range(1,13))
if metric=='SIV':
    mean_obs=np.mean(PIOMAS[(syear-1979):(eyear-1979)], axis=0)*1e12

years=[2312, 2435, 2413, 2356, 2387, 2332] #from low to high SIV anomaly #BAR = 2332, 2312, 2413, 2356, 2387, 2435 
ensembles=["%03d" % i for i in range(1,13)]
#ensembles=["001"]
run_length=1 #years
daily_matrix=np.zeros((6,12,1095))
days=365*3
smonth=["%02d" % i for i in range(1,12,2)]
emonth=["%02d" % i for i in range(0,12,2)]
emonth[0]="12"

monthly_matrix=np.zeros((6,(len(years)*12),36))


#Select observation data set
Data="Had2CIS"

#obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', str.lower(metric)).squeeze()

#Select Year Range
obs_years=range(1980,2019)
period=max(obs_years) - min(obs_years)

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
    obs2mask=obsin[((min(obs_years)-1980)*12):((max(obs_years)+2-1980)*12)]
    obs=np.reshape(obs2mask, ((period+2),12)).transpose()
else:
    obsingeo=np.delete(nc.getvar('/home/josmarti/Data/Observations/had2cis_128_64_195901_202004_sic.nc', 'SICN').squeeze(), 128, 2)
    gridpoint=nc.getvar('/home/josmarti/Data/areacella_fx_CanCM4_decadal2001_r0i0p0.nc','areacella')
    obs1=np.multiply(obsingeo,gridpoint)
    obsnh=np.delete(obs1, range(32), axis=1)
    obstemp=np.mean(np.mean(obsnh, axis=1), axis=1)
    obs2mask=obstemp[((min(obs_years)-1959)*12):((max(obs_years)+2-1959)*12)] 
    if Data == "Had2CIS":
        if metric == "SIA":
            obs=np.reshape(obs2mask, ((period+2),12)).transpose()
        elif metric == "SIE":
            obsin=nc.getvar('/home/josmarti/Data/Observations/had2cisSIE.nc', 'SICN').squeeze()[((min(obs_years)-1959)*12):((max(obs_years)+2-1959)*12)] 
            obs=np.reshape(obsin, ((period+2),12)).transpose()
    elif Data == "NSIDC":
        obsin=nc.getvar('/home/josmarti/Data/Observations/NSIDC_1979_2010_nh_siea.nc', str.lower(metric)).squeeze()
        obs=np.delete((np.reshape(obsin, ((period+3),12)).transpose()),range(min(obs_years)-1979),1)

original_obs=obs #used for standard deviation and end plot

for i in range(6):
    if i == 0:
        jd=0
    else:
        jd=1
    for e in ensembles:
        for year in years:
            if metric=='SIE':
                obs=annual_mean_obs
                if r!=0:
                    geosum1=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/%s/SIE_PP_run_full_%s_e%s_i%s%s_%s%s.nc' % (region,region,e,str(year),smonth[i],str(year+2+jd),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
                else:
                    geosum1=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PP_run_full_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+2+jd),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
                #geosum2=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PPrun_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+1),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
                #geosum3=nc.getvar('/home/josmarti/Data/PP_run_full/SIE/SIE_PPrun_e%s_i%s%s_%s%s.nc' % (e,str(year),smonth[i],str(year+2),emonth[i]), 'sicn').squeeze()*2*math.pi*6.371**2
            if metric=='SIA':
                obs=annual_mean_obs
                var1=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year),str(year))), 'sicn'), gridpoint)
                var2=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year+1),str(year+1))), 'sicn'), gridpoint)
                var3=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sicn.nc' % (e,str(year),str(year+2),str(year+2))), 'sicn'), gridpoint)
            if metric=='SIV':
                obs=np.mean(PIOMAS[(syear-1979):(eyear-1979)], axis=1)*1e12
                var1=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year),str(year))), 'sic'), gridpoint/913)
                var2=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+1),str(year+1))), 'sic'), gridpoint/913)
                var3=np.multiply(nc.getvar(('/home/josmarti/Data/PP_run_full/sc_dhfp1e_e%s_i%s_m01_%s_m01_%s_m12_sic.nc' % (e,str(year),str(year+2),str(year+2))), 'sic'), gridpoint/913)
            #var=nc.getvar('/home/josmarti/Data/PPrun_Jan/sc_dhfp1e_e%s_i%s_m01_daily_sicn.nc4' % (e,str(year)),'sicn')
            if metric != 'SIE':
                var1nh=np.delete(var1, range(32), axis=1)
                var2nh=np.delete(var2, range(32), axis=1)
                var3nh=np.delete(var3, range(32), axis=1)
                geosum1=np.sum(np.sum(var1nh, axis=1), axis=1)
                geosum2=np.sum(np.sum(var2nh, axis=1), axis=1)
                geosum3=np.sum(np.sum(var3nh, axis=1), axis=1)
                geosum1=np.concatenate((geosum1,geosum2,geosum3))
            #daily_matrix[years.index(year),ensembles.index(e),:]=np.sum(np.sum(var[:,32:64,:], axis=1), axis=1)
            if year==2413 and e=='011' and i==0: #i241501 e11 seems to have an extra timestep
                geosum1=geosum1[0:-1]
            monthly_matrix[i,(12*years.index(year)+ensembles.index(e)),:]=geosum1

#%%
"""for j in range(6): #MITCH DATA 1/3
    monthly_matrix[j,:,:]=np.reshape(np.transpose(loadmat('/home/josmarti/Downloads/SIE_16.mat')['metric_ensemble'][:,j,:,:], (0,2,1)), (72,36))
#"""
#%%

#clim_mean=data['control_clim'] #MITCH DATA 2/3
#control=np.reshape(data['metric_control'], (300,12)).transpose() #MITCH DATA 3/3
clim_mean=np.mean(control, axis=1, keepdims=True).transpose()
    
model=np.zeros((6,(len(years)*12),36))
truth=np.zeros((6,(len(years)*12),36))

ACC=np.zeros((12,36))
boot_temp=np.zeros((12,36), dtype=np.ndarray)
boot=np.zeros((12,36), dtype=np.ndarray)

if auto_corr:
    auto=np.zeros((12,36))
    boot_temp_auto=np.zeros((12,36), dtype=np.ndarray)
    boot_auto=np.zeros((12,36), dtype=np.ndarray)
    for i in range(12):
        for j in range(36):
            if (i-j-1)<0:
                auto[i,j]=np.corrcoef(control[(i-((j % 12)+1)),0:-((j/12)+1)],control[i,((j/12)+1)::])[1,0]
                boot_temp_auto[i,j]=bt.calc_corr_boot(control[(i-((j % 12)+1)),0:-((j/12)+1)],control[i,((j/12)+1)::], 1000)
                boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]
            else:
                if (j/12)==0:
                    auto[i,j]=np.corrcoef(control[(i-((j % 12)+1)),:],control[i,:])[1,0]
                    boot_temp_auto[i,j]=bt.calc_corr_boot(control[(i-((j % 12)+1)),:],control[i,:], 1000)
                    boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]
                else:
                    auto[i,j]=np.corrcoef(control[(i-((j % 12)+1)),0:-(j/12)],control[i,(j/12)::])[1,0]
                    boot_temp_auto[i,j]=bt.calc_corr_boot(control[(i-((j % 12)+1)),0:-(j/12)],control[i,(j/12)::], 1000)
                    boot_auto[i,j]=bt.calc_boot_stats(boot_temp_auto[i,j],sides=1,pval_threshold=0.05)[1]

for y in range(len(years)):
    for e in range(12):
        model[:,(y*12+e),:]=np.mean(np.delete(monthly_matrix[:,y*12:y*12+12,:], e, axis=1), axis=1)
        #truth[:,e,:]=monthly_matrix[]

truth=monthly_matrix

for i in range(0,12,2):
    for m in range(36):
        #model=np.repeat(np.mean(monthly_matrix[:,:,m], axis=1),12)
        #truth=np.reshape(monthly_matrix[:,:,m],72)
        if (m+i)<12:
            ACC[(i-12+m),m]=bt.calc_corr_choose(model[(i/2),:,m],truth[(i/2),:,m],clim_mean[0,(i-12+m)],clim_mean[0,(i-12+m)])
            boot_temp[(i-12+m),m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[(i-12+m),m]=bt.calc_boot_stats(boot_temp[(i-12+m),m],sides=1,pval_threshold=0.05)[1]
        elif (m+i)<24:
            ACC[(i-24+m),m]=bt.calc_corr_choose(model[(i/2),:,m],truth[(i/2),:,m],clim_mean[0,(i-24+m)],clim_mean[0,(i-24+m)])
            boot_temp[(i-24+m),m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[(i-24+m),m]=bt.calc_boot_stats(boot_temp[(i-24+m),m],sides=1,pval_threshold=0.05)[1]
        else:
            ACC[(i-36+m),m]=bt.calc_corr_choose(model[(i/2),:,m],truth[(i/2),:,m],clim_mean[0,(i-36+m)],clim_mean[0,(i-36+m)])
            boot_temp[(i-36+m),m]=bt.calc_corr_boot(model[(i/2),:,m],truth[(i/2),:,m],1000)
            boot[(i-36+m),m]=bt.calc_boot_stats(boot_temp[(i-36+m),m],sides=1,pval_threshold=0.05)[1]


#Interpolate ACCs
for j in range(1,12,2):
    ACC[j,0]=ACC[j-13,0]+0.5*(ACC[j-11,0]-ACC[j-13,0])
    boot_temp[j,0]=boot_temp[j-13,0]+0.5*(boot_temp[j-11,0]-boot_temp[j-13,0])
    boot[j,0]=bt.calc_boot_stats(boot_temp[j,0],sides=1,pval_threshold=0.05)[1]
    for n in range(2,35,2):
        ACC[j,n]=ACC[j,n-37]+0.5*(ACC[j,n-35]-ACC[j,n-37])
        boot_temp[j,n]=boot_temp[j,n-37]+0.5*(boot_temp[j,n-35]-boot_temp[j,n-37])
        boot[j,n]=bt.calc_boot_stats(boot_temp[j,n],sides=1,pval_threshold=0.05)[1]
    
for k in range(0,12,2):
    ACC[k,35]=ACC[k-1,35]+0.5*(ACC[k+1,35]-ACC[k-1,35])
    boot_temp[k,35]=boot_temp[k-1,35]+0.5*(boot_temp[k+1,35]-boot_temp[k-1,35])
    boot[k,35]=bt.calc_boot_stats(boot_temp[k,35],sides=1,pval_threshold=0.05)[1]
    for p in range(1,34,2):
        ACC[k,p]=ACC[k,p-1]+0.5*(ACC[k,p+1]-ACC[k,p-1])
        boot_temp[k,p]=boot_temp[k,p-1]+0.5*(boot_temp[k,p+1]-boot_temp[k,p-1])
        boot[k,p]=bt.calc_boot_stats(boot_temp[k,p],sides=1,pval_threshold=0.05)[1]
        

#%%         

#ACC=data['ACC_year_target']
jet3=colors.ListedColormap(loadmat('/home/josmarti/Downloads/cmap_jet3.mat')['cmap'], name='jet3') #Mitch's colors

if std_mask:
        if r!=0:
            if (len(np.std(original_obs, axis=1)[(np.std(original_obs, axis=1)*1e15)/np.sum(area) < 0.8])!=0 and std_mask):
                print ("STD mask in effect")
                for i in range(len(auto)):
                    if ((np.std(original_obs, axis=1)[i]*1e15)/np.sum(area)) < 0.8:
                        boot[i,:]=-1  # 0 will be seen as significant
                        ACC[i,:]=0


fig, ax = plt.subplots(figsize=(4,9))
if pstyle == "pc":
    #pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar') #CCCMA
    pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)  #Mitch's
    pc=rpl.add_pc(ax,range(13),range(37),ACC.transpose(),**pcparams)
elif pstyle == "cf":
    pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar',latlon=False)
    pc=rpl.add_cf(ax,range(1,13),range(1,37),ACC,**pcparams)
for month in range(len(boot[:,0])):
            for lead in range(len(boot[0,:])):
                if boot[month,lead]>0:   #MUST BE CHANGED BACK to >=!!!!!!
                    plt.scatter((month+0.5),(lead+0.5), color = 'black', s=20)
if thesis_figures:
    plt.title("%s PM" % region_titles[r], fontsize=16)
else:
    if r!=0:
        plt.title("%s Perfect Model Skill" % region, fontsize=16)
    else:
        plt.title("CanCM4 Perfect Model Skill", fontsize=16)
plt.colorbar(pc)
ax.invert_xaxis
ax.invert_yaxis
ax.set_xticks(np.arange(len(calendar.month_name[1:13]))+0.5)   
ax.set_xticklabels(calendar.month_name[1:13], rotation=90)    
ax.set_yticks(np.arange(36)+0.5)    
ax.set_yticklabels(range(36))  
plt.xlabel("Target Month")    
plt.ylabel("Lead (Months)")      
plt.show()            

if auto_corr:

    if std_mask:
        if r!=0:
            if (len(np.std(original_obs, axis=1)[(np.std(original_obs, axis=1)*1e15)/np.sum(area) < 0.8])!=0 and std_mask):
                print ("STD mask in effect")
                for i in range(len(auto)):
                    if ((np.std(original_obs, axis=1)[i]*1e15)/np.sum(area)) < 0.8:
                        boot_auto[i,:]=-1  # 0 will be seen as significant
                        auto[i,:]=0
    
    fig, ax = plt.subplots(figsize=(4,9))
    if pstyle == "pc":
        #pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar') #CCCMA
        pcparams=dict(clevs=np.arange(-0.8,1,0.1),cmap=jet3)  #Mitch's
        pc=rpl.add_pc(ax,range(13),range(37),auto.transpose(),**pcparams)
    elif pstyle == "cf":
        pcparams=dict(clevs=np.arange(-0.15,1.05,0.1),cmap='acccbar',latlon=False)
        pc=rpl.add_cf(ax,range(1,13),range(1,37),ACC,**pcparams)
    for month in range(len(boot[:,0])):
                for lead in range(len(boot[0,:])):
                    if boot_auto[month,lead]>0:   #MUST BE CHANGED BACK to >=!!!!!!
                        plt.scatter((month+0.5),(lead+0.5), color = 'black', s=20)
    if thesis_figures:
        plt.title("%s PM" % region_titles[r], fontsize=16)
    else:
        if r!=0:
            plt.title("%s PM Autocorrelation" % region, fontsize=16)
        else:
            plt.title("CanCM4 PM Autocorrelation", fontsize=16)
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

np.save('temp/%s.npy' % region, ACC.transpose()[0:12])
np.save('temp/%s_boot.npy' % region, boot_temp.transpose()[0:12])        
#%%

if show_timeseries:
    fig, ax = plt.subplots()
    for i in range(72): 
        plt.plot(range(12), monthly_matrix[0,i,0:12])
    if metric != "SIV":
        plt.plot(np.mean(original_obs*2*math.pi*6.371**2, axis=1), linewidth=2.5, color='k', label="Observations")
    else:
       plt.plot(mean_obs, linewidth=2.5, color='k', label="Observations") 
    plt.title("Montly mean SIV of first year run of 12 ensembles for six start years\n Observations (1979-2016)")
    plt.legend()
    ax.set_xticks(np.arange(len(calendar.month_name[1:13])))   
    ax.set_xticklabels(calendar.month_name[1:13], rotation=90, fontsize=12) 