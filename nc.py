# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 16:45:51 2015

@author: acrnrms
"""
import numpy as np # for array handling
from netCDF4 import Dataset

def openNC(filename):

    ncfile = Dataset(filename,'r')
    return ncfile
    
    
def getvar(filename,field,remcycl=False):
    ncfile = openNC(filename)
    ncvar = ncfile.variables[field][...]


    # Apply any scaling and offsetting needed:
    try:
        var_offset = ncvar.add_offset
        print 'var_offset ' + str(var_offset)

    except:
        var_offset = 0
            
    try:
        var_scale = ncvar.scale_factor
        print 'var_scale ' + str(var_scale)
    except:
        var_scale = 1

    # Read missing value
#    try:
#        miss=ncvar.missing_value
#        print 'missing_value' + str(miss)
#    except:
#        miss=-9e99
#        print 'no missing value'          
    fld = ncvar*var_scale + var_offset
    
#    fld = np.ma.masked_where(fld==miss,fld)

    # Remove last cyclic longitude point:
    if remcycl!=False:
       # remove extra lon
       fld = np.squeeze(fld[...,0:-1]) # @@



    ncfile.close()
    return fld
