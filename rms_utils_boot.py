# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:15:11 2016

@author: acrnrms

"""
import numpy as np
import matplotlib.pyplot as plt
testcode=False

##################################################################################
def calc_corr(x,y):
    """ Compute Correlation between x and y

    Parameters
    ----------
    x : timeseries 1 [time,...]
    y : timeseries 2 [time,...]

    Returns
    -------
    corr: correlation between x and y [can be multidimensional]

Created on Thu Oct 27 11:15:11 2016

@author: acrnrms
"""
    nt=np.shape(x)[0]
    #calculate anomalies
    x_a = x - x.mean(axis=0)
    y_a = y - y.mean(axis=0)
    #calculate covariance
    cov = np.sum(x_a*y_a,axis=0)/(nt-1)
    #calculate standard deviations
    std_x = np.sqrt(np.sum(x_a**2.,axis=0)/(nt-1))
    std_y = np.sqrt(np.sum(y_a**2.,axis=0)/(nt-1))
    #calculate rho (correlation coefficient)
    corr = cov/(std_x*std_y)
    #print corr
    if np.size(np.shape(x))==1:    
        if np.isnan(corr): corr=0.0
    else:
        corr[np.isnan(corr)==True] = 0.0
    return corr

##################################################################################
def calc_corr_boot(x, y,nboot):    
    """ Compute Bootstrapped Correlation between x and y

    Parameters
    ----------
    x : timeseries 1 [time,...]
    y : timeseries 2 [time,...]

    Returns
    -------
    corr: bootstrapped correlation between x and y [boot,...]
    """ 

    nt=np.shape(x)[0]  
    nboot=np.array(nboot,dtype=int)
    dims=np.concatenate(([nboot],np.array(np.shape(x)[1::])))
    #print (dims)
    #dims=np.int(dims)
    corr_boot=np.zeros(int(dims))
    for i in range(nboot):
        resample_i = np.floor(np.random.rand(nt)*nt).astype(int)
        #print resample_i
        corr_boot[i,...]=calc_corr(x[resample_i,...],y[resample_i,...])
    return corr_boot    

##################################################################################
def calc_mean_boot(x, nboot):    
    """Parameters
    ----------
    x : sample 1 
    """
    nx=np.shape(x)[0]
    mean_boot=np.zeros((nboot,1))   
    for i in range(nboot):
        resample_i = np.floor(np.random.rand(nx)*nx).astype(int)
        mean_boot[i]=np.mean(x[resample_i])

    return mean_boot     
     


##################################################################################
def calc_boot_stats(x_boot,sides=1,pval_threshold=0.05):
    """ Calculate statistics from a bootstrap

    Parameters
    ----------
    x_boot : Bootstrap [nboot,...]

    Returns
    -------
    x_boot_mean:  Mean [...]
    x_boot_min: 5th percentile (or 2.5 if sides=2)
    x_boot_max: 95th percentile (or 97.5 if sides=2)    
    x_pval: pval that x >0     
    
    """ 
    #dims
    nboot=np.shape(x_boot)[0]      
    #pval_threshold_pencentile
    if sides==1: pval_threshold_pctl=pval_threshold*100
    if sides==2: pval_threshold_pctl=pval_threshold*100/2
    #mean
    x_boot_mean=np.mean(x_boot,axis=0)
    #min/max
    x_boot_min=np.percentile(x_boot,pval_threshold_pctl,axis=0)
    x_boot_max=np.percentile(x_boot,100-pval_threshold_pctl,axis=0)
    
    #pval
    x_pval=(np.sort(x_boot,axis=0)>0).argmax(axis=0)/(np.float(nboot)*np.float(sides))

    return x_boot_mean,x_boot_min,x_boot_max,x_pval

##################################################################################

#####TEST################
if testcode:
    #random data
    nboot=10000
    x=np.random.rand(32)
    y=np.random.rand(32)
    x=np.zeros((32,1))
    y=np.zeros((32,1))
    
    #corr_boot
    corr_boot=calc_corr_boot(x,y,nboot)
    plt.hist(corr_boot)
    #corr_boot_stats
    corr=calc_corr(x,y) 
    print 'CORR: ' + str(corr)
    
    boot_mean=calc_boot_stats(corr_boot)[0]
    boot_min=calc_boot_stats(corr_boot)[1]
    pval=calc_boot_stats(corr_boot)[3]
    
    print 'boot_mean: '     + str(boot_mean)
    print 'boot_min: '      + str(boot_min)

    print 'pval: '          + str(pval) 
    print 'CONF INTERVAL: ' + str(boot_mean-boot_min)
    
##############################################################################
    
def calc_corr_choose(x,y,x_mean,y_mean):
    """ Compute Correlation between x and y

    Parameters
    ----------
    x : timeseries 1 [time,...]
    y : timeseries 2 [time,...]
    x_mean:
    y_mean:

    Returns
    -------
    corr: correlation between x and y [can be multidimensional]

Created on Thu Oct 27 11:15:11 2016

@author: acrnrms
"""
    nt=np.shape(x)[0]
    #calculate anomalies
    x_a = x - x_mean
    y_a = y - y_mean
    #calculate covariance
    cov = np.sum(x_a*y_a,axis=0)/(nt-1)
    #calculate standard deviations
    std_x = np.sqrt(np.sum(x_a**2.,axis=0)/(nt-1))
    std_y = np.sqrt(np.sum(y_a**2.,axis=0)/(nt-1))
    #calculate rho (correlation coefficient)
    corr = cov/(std_x*std_y)
    #print corr
    if np.size(np.shape(x))==1:    
        if np.isnan(corr): corr=0.0
    else:
        corr[np.isnan(corr)==True] = 0.0
    return corr
    
##################################################################################
def calc_corr_boot_choose(x, y, mean, nboot):    
    """ Compute Bootstrapped Correlation between x and y

    Parameters
    ----------
    x : timeseries 1 [time,...]
    y : timeseries 2 [time,...]

    Returns
    -------
    corr: bootstrapped correlation between x and y [boot,...]
    """ 

    nt=np.shape(x)[0]  
    nboot=np.array(nboot,dtype=int)
    dims=np.concatenate(([nboot],np.array(np.shape(x)[1::])))
    #print (dims)
    #dims=np.int(dims)
    corr_boot=np.zeros(int(dims))
    for i in range(nboot):
        resample_i = np.floor(np.random.rand(nt)*nt).astype(int)
        #print resample_i
        corr_boot[i,...]=calc_corr_choose(x[resample_i,...],y[resample_i,...], mean, mean)
    return corr_boot    



