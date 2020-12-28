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
import cdo#; cdo = cdo.Cdo()
import matplotlib.pyplot as plt # for basic plotting
import rms_plots as rpl
import nc as  nc

###############################################################################
#######READ########## #########################################################
###############################################################################
infile = '/home/josmarti/Data/iceregions_128_64.nc'
lon=nc.getvar(infile,'lon')
lat=nc.getvar(infile,'lat')
region=nc.getvar(infile,'REG').squeeze()
regionlabs=['-1ocean','0land','1ARC','2BER','3STL','4BAF','5GRE','6BAR','7KAR','8LAP','9ESI','10CHU','11BEA','12CAN','13HUD','14OKH']
###############################################################################
#######PLOT####################################################################
###############################################################################
###colors####################################################################
clevs=np.arange(17)-1.5
clevslab=np.arange(17)-1

full_cmap = plt.get_cmap('jet')
colors=full_cmap([0,20,40,60,80,100,110,120,130,140,160,180,200,220,240,255])

####open figs####################################################################
ncol=1; nrow=1; fig, axs = plt.subplots(nrow,ncol); fig.set_size_inches(ncol*2,nrow*2.75)
fig.subplots_adjust(left=0.1,right=0.81)
fig.suptitle('Navy subregions' ,fontsize=10,y=0.93)

####plot####################################################################
bm,pc=rpl.plot_latlon(lat,lon,region,axis=axs,maptype='nsnhz',shadetype='pcolor',clevs=clevs,colors=colors,extend='neither',supp_cb=True)
####cbar####################################################################
cax = fig.add_axes([0.82,0.13,0.03,0.75])
cbar=fig.colorbar(pc,cax=cax,orientation='vertical',ticks=clevslab)
cbar.ax.set_yticklabels(regionlabs)
cbar.ax.tick_params(labelsize=6,pad=1,length=0.1)
######save####################################################################
#fig.savefig('regions.png',dpi=400)



