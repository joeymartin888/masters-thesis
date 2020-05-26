#!/usr/bin/env python
# coding: utf-8

# # CanESM5 GMD paper Figure 14x
# ## Northern and Southern Annular modes, based on psl EOF 1.
# 
# CMOR variables: `psl`.
# 
# ***This version has maps of the mean from CanESM5 and anomalies from observations***
# 
# 
# ### history
# - NCS, setup basic example, 2019-05-08

# In[229]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
get_ipython().magic(u'matplotlib inline')
import xarray as xr
import numpy as np
import os, sys
from pprint import pprint
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy 
from datetime import datetime
import cmocean
from cdo import Cdo
from eofs.xarray import Eof
import os
os.environ["PATH"] += ':/home/ords/crd/ccrn/scrd104/miniconda3/bin/'
cdo = Cdo()
cdo = Cdo()
cdo.setCdo('/home/ords/crd/ccrn/scrd104/miniconda3/bin/cdo')
import cartopy
import matplotlib.path as mpath


# In[230]:


runCDO = False

if runCDO:   
    CanESM5_psl = 'psl_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc'
    ifile = os.path.join('/space/hall2/sitestore/eccc/crd/CMIP6/final/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p1f1/Amon/psl/gn/v20190306/', CanESM5_psl)
    ofile = 'processed_data/remap-woa09_' + CanESM5_psl
    cdo.remapbil('obs/woa/woa09/uncs_woa09_ann_tpot.nc', input=ifile, output=ofile)
    
    ERAInt_psl = 'psl_Amon_ERA-Int.nc'
    ifile = os.path.join('obs/ERA-Int/', ERAInt_psl)
    ofile = 'processed_data/remap-woa09_' + ERAInt_psl
    cdo.remapbil('obs/woa/woa09/uncs_woa09_ann_tpot.nc', input=ifile, output=ofile)


# In[245]:


start = '1981-01-01'
end = '2014-12-31'
startdate = datetime.strptime(start,'%Y-%m-%d')
enddate   = datetime.strptime(end,'%Y-%m-%d')

psl = xr.open_dataset('processed_data/remap-woa09_psl_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc')['psl'] #units Pa
psl = psl.sel(time=slice(start, end))

psl_obs  = xr.open_dataset('processed_data/remap-woa09_psl_Amon_ERA-Int.nc')['psl'] #units Pa
psl_obs  = psl_obs.sel(time=slice(start, end))


# In[246]:


psl_sof20s = psl.sel(lat=slice(-90,-20))
psl_sof20s =psl_sof20s - psl_sof20s.mean(dim='time')
coslat = np.cos(np.deg2rad(psl_sof20s.coords['lat'].values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
#psl_sof20s
solver = Eof(psl_sof20s, weights=wgts)
sh_eof = solver.eofsAsCorrelation(neofs=1)
var_s = solver.varianceFraction(neigs=1)

psl_sof20s_obs = psl_obs.sel(lat=slice(-90,-20))
psl_sof20s_obs = psl_sof20s_obs - psl_sof20s_obs.mean(dim='time')
#psl_sof20s
solver_obs = Eof(psl_sof20s_obs, weights=wgts)
sh_eof_obs = solver_obs.eofsAsCorrelation(neofs=1)
var_s_obs = solver_obs.varianceFraction(neigs=1)


# In[247]:


import iris
import iris.coord_categorisation

cube = iris.load_cube('processed_data/remap-woa09_psl_Amon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc')
cube2 = iris.load_cube('processed_data/remap-woa09_psl_Amon_ERA-Int.nc')

iris.coord_categorisation.add_season(cube, 'time', name='clim_season')
iris.coord_categorisation.add_season_year(cube, 'time', name='season_year')

iris.coord_categorisation.add_season(cube2, 'time', name='clim_season')
iris.coord_categorisation.add_season_year(cube2, 'time', name='season_year')


annual_seasonal_mean = cube.aggregated_by(['clim_season', 'season_year'],iris.analysis.MEAN)
annual_seasonal_mean2 = cube2.aggregated_by(['clim_season', 'season_year'],iris.analysis.MEAN)


# In[248]:


psl = xr.DataArray.from_iris(annual_seasonal_mean.extract(iris.Constraint(clim_season='djf'))).sel(time=slice(start, end))
psl_obs = xr.DataArray.from_iris(annual_seasonal_mean2.extract(iris.Constraint(clim_season='djf'))).sel(time=slice(start, end))


# In[249]:


psl_nof20n = psl.sel(lat=slice(20,90))
psl_nof20n = psl_nof20n - psl_nof20n.mean(dim='time')

coslat = np.cos(np.deg2rad(psl_nof20n.coords['lat'].values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
#psl_sof20s
nsolver = Eof(psl_nof20n, weights=wgts)
nh_eof = nsolver.eofsAsCorrelation(neofs=1)
var_n = nsolver.varianceFraction(neigs=1)

psl_nof20n_obs = psl_obs.sel(lat=slice(20,90))
psl_nof20n_obs = psl_nof20n_obs - psl_nof20n_obs.mean(dim='time')

#psl_sof20s
nsolver_obs = Eof(psl_nof20n_obs, weights=wgts)
nh_eof_obs = nsolver_obs.eofsAsCorrelation(neofs=1)
var_n_obs = nsolver_obs.varianceFraction(neigs=1)


# In[259]:


fig=plt.figure(figsize=(12,8))
   
# CanESM5
ax = plt.subplot2grid((2, 2), (0, 0),projection=ccrs.NorthPolarStereo())
ax.pcolormesh(psl_nof20n.lon, psl_nof20n.lat,nh_eof[0,...],vmin=-1,vmax=1,cmap=plt.cm.RdBu_r,transform=ccrs.PlateCarree())
ax.contour(psl_nof20n.lon, psl_nof20n.lat,nh_eof[0,...],10,transform=ccrs.PlateCarree(),colors='k', zorder=20, alpha=0.4)

ax.set_extent([-180, 180, 20, 90],ccrs.PlateCarree())
ax.coastlines()
#ax.gridlines()
#ax.set_boundary(circle, transform=ax.transAxes)
#ax.add_feature(cartopy.feature.LAND, zorder=10, edgecolor='grey')
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_title('(a) CanESM5 NAM ({0:.2f})'.format(var_n.values[0]))

# Had
ax = plt.subplot2grid((2, 2), (0, 1),projection=ccrs.NorthPolarStereo())
ax.pcolormesh(psl_nof20n_obs.lon, psl_nof20n_obs.lat,-nh_eof_obs[0,...],vmin=-1,vmax=1,cmap=plt.cm.RdBu_r,transform=ccrs.PlateCarree())
ax.contour(psl_nof20n.lon, psl_nof20n.lat,-nh_eof_obs[0,...],10,transform=ccrs.PlateCarree(),colors='k',zorder=20, alpha=0.4)
ax.set_extent([-180, 180, 20, 90],ccrs.PlateCarree())
ax.coastlines(linewidth=1)
#ax.gridlines()
#ax.set_boundary(circle, transform=ax.transAxes)
#ax.add_feature(cartopy.feature.LAND, zorder=10, edgecolor='grey')
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_title('(b) ERA-Int NAM ({0:.2f})'.format(var_n_obs.values[0]))

# CanESM5
ax = plt.subplot2grid((2, 2), (1, 0),projection=ccrs.SouthPolarStereo())
ax.pcolormesh(psl_sof20s.lon, psl_sof20s.lat,-sh_eof[0,...],vmin=-1,vmax=1,cmap=plt.cm.RdBu_r,transform=ccrs.PlateCarree())
ax.contour(psl_sof20s.lon, psl_sof20s.lat,-sh_eof[0,...],transform=ccrs.PlateCarree(),colors='k',zorder=20, alpha=0.4)

ax.set_extent([-180, 180, -90, -20],ccrs.PlateCarree())
ax.coastlines()
#ax.gridlines()
#ax.set_boundary(circle, transform=ax.transAxes)
#ax.add_feature(cartopy.feature.LAND, zorder=10, edgecolor='grey')
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_title('(c) CanESM5 SAM ({0:.2f})'.format(var_s.values[0]))

# Had
ax = plt.subplot2grid((2, 2), (1, 1),projection=ccrs.SouthPolarStereo())
ax.pcolormesh(psl_sof20s_obs.lon, psl_sof20s_obs.lat,-sh_eof_obs[0,...],vmin=-1,vmax=1,cmap=plt.cm.RdBu_r,transform=ccrs.PlateCarree())
ax.contour(psl_sof20s_obs.lon, psl_sof20s_obs.lat,-sh_eof_obs[0,...],transform=ccrs.PlateCarree(),colors='k',zorder=20, alpha=0.4)

ax.set_extent([-180, 180, -90, -20],ccrs.PlateCarree())
ax.coastlines()
#ax.gridlines()
#ax.set_boundary(circle, transform=ax.transAxes)
#ax.add_feature(cartopy.feature.LAND, zorder=10, edgecolor='grey')
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.set_title('(d) ERA-Int SAM ({0:.2f})'.format(var_s_obs.values[0]))
plt.subplots_adjust(wspace=-0.5)
fig.savefig('plots/fig-14x_atmos_nam_sam.png', bbox_inches='tight')


# In[251]:


np.corrcoef(nh_eof[0,...].values.flatten(), nh_eof_obs[0,...].values.flatten())


# In[252]:


np.corrcoef(sh_eof[0,...].values.flatten(), sh_eof_obs[0,...].values.flatten())


# In[ ]:




