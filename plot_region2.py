#####################################################
##
## PROGRAM: calc_fudifd_acc_clim_region.sh
##
## 21 Oct 2015: M. Sigmond 
##
## PURPOSE: Calculate Regional means of climatologies and skill FUD and IFD
##
## AUTHOR: Michael Sigmond (CCCma)
##
## INPUT 
#####################################################
import numpy as np
import cdo; cdo = cdo.Cdo()
import matplotlib.pyplot as plt # for basic plotting
import rms_plots as rpl
import nc as  nc

###############################################################################
#######READ########## #########################################################
###############################################################################
infile = '/HOME/rms/DATA/MASKS/ICE/NAVY/iceregions_128_64.nc'
lon=nc.getvar(infile,'lon'); nlon=len(lon)
lat=nc.getvar(infile,'lat'); nlat=len(lat)
region=nc.getvar(infile,'REG').squeeze()
regionlabs=['-1ocean','0land','1ARC','2BER','3STL','4BAF','5GRE','6BAR','7KAR','8LAP','9ESI','10CHU','11BEA','12CAN','13HUD','14OKH']
###############################################################################
#######PLOT####################################################################
###############################################################################
###colors####################################################################
clevs=np.arange(17)-1.5
clevslab=np.arange(17)-1


####open figs####################################################################
ncol=1; nrow=1; fig, axs = plt.subplots(nrow,ncol); fig.set_size_inches(ncol*2,nrow*2.75)
fig.subplots_adjust(left=0.1,right=0.81)
#fig.suptitle('Navy subregions' ,fontsize=10,y=0.93)

####plot####################################################################
bm=rpl.make_basemap(maptype='nsnhz',axis=axs) 
full_cmap = plt.get_cmap('bluegrayred11')
icols=[7,7,2,7,2,7,7,2,7,2,2]
for ii,ireg in enumerate([2,4,5,6,7,8,9,10,11,13,14]):
    print icols[ii]
    colors=full_cmap([0,icols[ii]]);clevs=[0.5]
    y=np.zeros([nlat,nlon]);
    y[region==ireg]=1    
    pc=rpl.add_plot_latlon(lat,lon,np.ma.masked_where(y<1,y),bm,axis=axs,clevs=clevs,colors=colors,supp_cb=True,shadetype='pcolor') 

axs.text(0.13,0.24,'Hudson',fontsize=5,color='k',transform = axs.transAxes)
#axs.text(0.15,0.21,'Bay',fontsize=5,transform = axs.transAxes)
axs.text(0.35,0.24,'Baffin',fontsize=5,color='k',transform = axs.transAxes)
#axs.text(0.35,0.21,'Bay',fontsize=5,transform = axs.transAxes)
axs.text(0.32,0.12,'Labrador',fontsize=5,color='k',transform = axs.transAxes)
#axs.text(0.35,0.11,'Sea',fontsize=5,transform = axs.transAxes)
axs.text(0.55,0.32,'Greenland',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.67,0.45,'Barents',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.65,0.54,'Kara',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.49,0.60,'Laptev',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.35,0.64,'E.Siberian',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.22,0.60,'Chukchi',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.18,0.50,'Beaufort',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.10,0.70,'Bering',fontsize=5,color='k',transform = axs.transAxes)
axs.text(0.33,0.89,'Okhotsk',fontsize=5,color='k',transform = axs.transAxes)


####cbar####################################################################

######save####################################################################
fig.savefig('regions2.png',dpi=400)
fig.savefig('regions2.pdf',dpi=400)


