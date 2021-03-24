#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 23:12:20 2021

@author: josmarti
"""

import numpy as np

def Mitch_detrend(data):
    
    ############
    #
    #
    #
    ############3
    
    x=np.arange(data.shape[-1])
    ret=np.zeros(data.shape,)
    
    if len(ret.shape) == 2:
        for t in range(ret.shape[0]):
            for year in range((len(x)-3)):
                past=data[t,0:(year+3)]
                A = np.vstack([x[0:(year+3)], np.ones(len(past))]).T
                m, c = np.linalg.lstsq(A, data[t,0:(year+3)], rcond=None)[0]
                ret[t,(year+3)]=data[t,(year+3)]-(m*x[year+3]+c)
            for year in range(2):
                ret[t,year]=data[t,year]-np.mean(data[t,0:year+1])
    
    if len(ret.shape) == 3:
        for t in range(ret.shape[0]):
            for l in range(ret.shape[1]):
                for year in range((len(x)-3)):
                    past=data[l,t,0:(year+3)]
                    A = np.vstack([x[0:(year+3)], np.ones(len(past))]).T
                    m, c = np.linalg.lstsq(A, data[l,t,0:(year+3)], rcond=None)[0]
                    ret[l,t,(year+3)]=data[l,t,(year+3)]-(m*x[year+3]+c)
                for year in range(2):
                    ret[l,t,year]=data[l,t,year]-np.mean(data[l,t,0:year+1])
    
    return ret
    
        
    
    