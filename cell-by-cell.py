#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:06:07 2021

@author: josmarti
"""
import netCDF4 as nc
import numpy as np
import calendar

"""s = pxssh.pxssh()
if not s.login ('aeolus.seos.uvic.ca', 'josmarti', open("/home/josmarti/pword.txt").read(), port=11100):
    print "SSH session failed on login."
    print str(s)
else:
    print "SSH session login successful"
    s.sendline ('ls -l')
    s.prompt()         # match the prompt
    print s.before     # print everything before the prompt.
    s.logout()"""

data=np.zeros(((12*(2018-1980)),12,180,360))
months=["%02d" % i for i in range(1,13)]

for i in range((12*(2018-1980))):
     data[i,:,:,:]=np.mean(nc.getvar(('/home/josmarti/Data/NEW/sic_monthly_CCCma-CanCM4_NEW_1x1_grid_i%i%s01.nc' % ((1980+(i/12)),months[(i % 12)])), 'sic').squeeze(), axis=1)


#%%

"""test=np.zeros((456,12)).astype(str)

for m in range(12):
    for n in range(m,456,12):
        for i in range(12):
            test[n,i]="%s%i" % (calendar.month_name[1:13][(m % 12)+i-12],i)"""

new_frame=np.zeros((12,12,38,180,360)).astype(str)

for target in range(12):
    for lead in range(12):
        for i in range(target-lead,456,12):
            new_frame[lead,target]=data[i,lead,:,:]

np.save('new_frame',new_frame)