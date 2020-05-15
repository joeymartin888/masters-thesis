#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:02:45 2020

@author: josmarti
"""
import numpy as np
#from mpl_toolkits.basemap import Basemap #(not installed)
import matplotlib.pyplot as plt # for basic plotting
import mpl_toolkits as mpltk
import matplotlib.colors as col
from matplotlib.colors import from_levels_and_colors
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

def centre_to_edge(lonin,latin):
    ###LON
    nlon=len(lonin);lonout=np.zeros(nlon);
    lonout[0]=(lonin[-2]+360+lonin[0])/2  
    for ilon in range(1,nlon):
        lonout[ilon]=(lonin[ilon]+lonin[ilon-1])/2
    ###LAT
    nlat=len(latin); latout=np.zeros(nlat+1);
    if latin[1]>latin[0]: #lat increasing [S to N]
        latout[0]=latin[0]-(latin[1]-latin[0])/2.    
        for ilat in range(1,nlat):
            latout[ilat]=(latin[ilat]+latin[ilat-1])/2  
        latout[nlat]=89.99
    else:
        latout[0]=89.99   
        for ilat in range(1,nlat):
            latout[ilat]=(latin[ilat]+latin[ilat-1])/2  
        latout[nlat]=latin[-1]-(latin[-2]-latin[-1])/2.    
        
    return lonout,latout  

def make_cartopy(projection=ccrs.LambertCylindrical(), figsize=(6, 4), resolution='110m'):
    ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection=projection))
    ax.set_global()
    ax.coastlines(resolution=resolution, color='k')
    # Only PlateCarree and Mercator plots are currently supported.
    gl = ax.gridlines(draw_labels=False)
    ax.add_feature(cfeature.LAND, facecolor='0.75')
    return ax

def add_plot_latlon_cartopy(lat,lon,fld,map,axis=None,cint='',clevs='',colors='',
                shadetype='contourf',supp_cb='False',extend='both'): 

    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: add_colorplot_latlon
    """      
    ##### Build grid, add cyclic lon if needed, shift lons+lats in case of pcolor plot
    if np.mod(lon.shape,2) == 0:
       fld,lon = cutil.add_cyclic_point(fld,lon) #changed to cartopy by Joey Martin
    if shadetype == 'pcolor':
        lon_edge,lat_edge=centre_to_edge(lon,lat)
        lons, lats = np.meshgrid(lon_edge,lat_edge) ###NB: len(lat_edge)=nlat+1
        fld = np.ma.masked_invalid(fld) #http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values                
    else:
        lons, lats = np.meshgrid(lon,lat)  
    #############Colors
    if colors=='':
        full_cmap = plt.get_cmap('bluegrayred19')
        colors=full_cmap(np.arange(19))

    #############clevs###################################### 
    # options:
    # 1) give clevs --> clevs=clevs
    # 2) give cint --> clevs are constructed centred around 0
    # 3) no clevs --> 19 autolevels

    ncolors=colors.shape[0]; #nr of colors in cbar
    if cint !='':  #construct from cint (if given)
       crange=cint*(ncolors-2)/2; cmin=-crange; cmax=crange
       clevs=np.linspace(cmin,cmax,ncolors-1);
    #############CMAP######################################
    if clevs != '':        
       cmaploc, norm = from_levels_and_colors(clevs, colors ,extend=extend)
       pcparams=dict(latlon=True,cmap=cmaploc,norm=norm)
    else:
       cmaploc=col.ListedColormap(colors)         
       pcparams=dict(latlon=True,cmap=cmaploc)
    #############plot######################    
    if shadetype == 'pcolor':
#        print pcpararms
        pc=plt.pcolormesh(lons,lats,fld,**pcparams)
    elif shadetype == 'contourf':
         if clevs != '':
             pcparams['levels']=clevs
             pcparams['extend']='both'
             pc=map.contourf(lons,lats,fld,**pcparams)
         else:
             pc=map.contourf(lons,lats,fld,19,**pcparams)    
    else:
        print "Incorrect shadetype. Choose pcolor or contourf"
        return -1
    ############cbar######################
    if supp_cb==False:
        if clevs != '':
            map.colorbar(pc,location='right',pad="5%",ticks=clevs)
        else:
            map.colorbar(pc,location='right',pad="5%")
    return pc


#_ = ax.plot(lon, lat, transform=ccrs.Geodetic(), **kw)