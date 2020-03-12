"""
Created on Wed Sep 16 18:21:18 2015

@author: acrnrms
"""
import numpy as np
#from mpl_toolkits.basemap import Basemap (not installed)
import matplotlib.pyplot as plt # for basic plotting
import mpl_toolkits as mpltk
import matplotlib.colors as col
from matplotlib.colors import from_levels_and_colors
import matplotlib.cm as cm

############################################################################
#COLOR MAPS#################################################################
############################################################################
def register_rms_cmaps(cmap='all'):
    """create my personal colormaps with discrete colors and register them.
       default is to register all of them. can also specify which one.
    User defined:
        bluegrayred19: default
        blue0red19: 
        blue0red11giss: Mimicing GISS T-trends (used for hiatus studies)
        bluegrayred9dark:  Used for retreat/advance skill studies
        bluegrayred11dark: Used for retreat/advance skill studies 
        
        


    """
    print 'registering cmaps'

#bluegrayred19#########
    # blueish at top, gray in middle, reddish at bottom
    colors = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [200, 200, 200],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)/255
    
    thecmap = col.ListedColormap(colors,'bluegrayred19')
    cm.register_cmap(cmap=thecmap)

#blue0red19#########
    # As bluegrayred19, but white in middle

    colors = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [255, 255, 255],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)/255
    
    thecmap = col.ListedColormap(colors,'blue0red19')
    cm.register_cmap(cmap=thecmap)

#bluegray0red20#########
    # As bluegray0red20, but white added

    colors = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [255, 255, 255],\
                       [225, 225, 225],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)/255
    
    thecmap = col.ListedColormap(colors,'blue0grayred20')
    cm.register_cmap(cmap=thecmap)



#blue0red11giss#########
    # blueish at top, white in middle, yellow and red at bottom
    # Mimicing GISS temp colors, but yellow less looking like pee as Fyfe's request
    colors = np.array([ [131,63,233], \
                       [71,137,252], \
                       [125,206,253],\
                       [165,250,255],\
                       [213,255,226],\
                       [255,255,255],\
                       [255,255,200],\
                       [255,210,27],\
                       [250,173,19],\
                       [255,0,0],\
                       [132,30,30]],\
                       dtype=float)/255.

    thecmap = col.ListedColormap(colors,'blue0red11giss')
    cm.register_cmap(cmap=thecmap)


#bluegrayred11dark#########
    # blueish at top, gray in middle, reddish at bottom
    # Adapted 11-class RdBu from colorbrewer2.org:
    # Recipe:
    # 1) 11-class Rdbu
    # 2) Replace the white color with gray (colorblind friendly)
    colors = np.array([ [5,48,97], \
                       [33,102,172], \
                       [67,147,195],\
                       [146,197,222],\
                       [209,229,240],\
                       [130,130,130],\
                       [253,219,199],\
                       [244,165,130],\
                       [214,96,77],\
                       [178,24,43],\
                       [103,0,31]],\
                    dtype=float)/255.

    thecmap = col.ListedColormap(colors,'bluegrayred11dark')
    cm.register_cmap(cmap=thecmap)

#bluegrayred9dark#########
    # Adapted 11-class RdBu from colorbrewer2.org:
    # Recipe:
    # 1) Pick the 4 darkest and 4 lightest colors from 11-class Rdbu
    # 2) Replace the 3 middle ones with a gray shading (colorblind friendly)
    colors = np.array([ [5,48,97], \
                       [33,102,172], \
                       [67,147,195],\
                       [146,197,222],\
                       [130,130,130],\
                       [244,165,130],\
                       [214,96,77],\
                       [178,24,43],\
                       [103,0,31]],\
                    dtype=float)/255.

    thecmap = col.ListedColormap(colors,'bluegrayred9dark')
    cm.register_cmap(cmap=thecmap)

#sic#########
    colors = np.array([[9, 60, 112],\
                      [255, 255, 255]],\
                      dtype=float)/255.
             
    thecmap = col.ListedColormap(colors,'sic')
    cm.register_cmap(cmap=thecmap)




#
register_rms_cmaps()

############################################################################
#AXIS FOR FIGURES###########################################################
############################################################################

###################################################
def make_bm(ax,
             region='glob',latnps0=50,latsps0=-50.,lonps0=270.,
             coastlinewidth=0.5,coastlinecolor='black',
             landmask=False,oceanmask=False,                   
             drawgrid=False,
             landfillcol=[0.85,0.85,0.85],
             oceanfillcol=np.array([176, 237, 245])/255.):

    """Make a basemap for latlon plot. 


    Parameters:
    -----------
      ax : [axis]
            Subplot axis

      *region* : ['glob' | string] 
                Region for the plots. Current options:
                #GLOBAL:
                'glob': Global, standard (cylindrical) projections 
                'glob_rob': Global, Robinson projection
                #POLAR:
                'nps': North Pole stereographic (latnps0,lonps0)
                'sps': South Pole stereographic (latsps0,lonps0)
                'nps2': As nps, but rectagle like NSIDC
                'nps3': As nps2, but zoomed in
                'baf' : zoom on Baffin
                #SUBREGIONS:
                'nh': NH, standard (cylindrical) projections 
                'tpo_na': Tropical Pacific Ocean + North America (cylindrical) 
                'atl': Atlantic (cylindrical) 
                'na': North Atlantic (cylindrical) 

      *latnps0* : [ *50.* | float ]         
                Latitude boundary for nps plots
                
      *latsps0* : [ *-50.* | float ]         
                Latitude boundary for nps plots

                
      *lonps0* : [ *270.* | float ]         
                Longitude centre for nps and sps plots

      *coastlinewidth* : [ *0.5* | float ]         
                Coastline width

      *coastlinecolor* : [ *'black'* | string ]         
                Coastline color

      *landmask*: [ *'False'* | Boolean ]
                If True: Fill continents and lakes with gray color  

      *oceanmask*: [ *'False'* | Boolean ]
                If True: Drawmapboundary with coral color
      *drawgrid*: [ *'False'* | Boolean ]
                If True:Draw grid every 30 degrees 

      *landfillcol*=[ *[0.85,0.85,0.85]* | arry ]
                Fill color of land if landmask=True

      *oceanfillcol*=[ np.array([176, 237, 245])/255. | array ]
                Fill color of land if oceanmask=True

      Returns:
      --------
          bm: Basemap handle
           
    """      
    polythres=50
    ##########global
    if region == 'glob': #latlon  (cylindrical)         
         lat1=-90; lat2=90;lon1=-180; lon2=180
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif region == 'glob_rob':
       mapparams = dict(projection='robin',lon_0=0,resolution='l')
    ##########Poles
    elif region == 'nps': #North Pole stereographic 
         mapparams = dict(projection='npstere',boundinglat=latnps0,lon_0=lonps0,
                     resolution='c',round=True)
         if latnps0==50: polythres=35
         if latnps0==40: polythres=55            
    elif region == 'sps': #South Pole stereographic 
         mapparams = dict(projection='spstere',boundinglat=latsps0,lon_0=lonps0,
                     resolution='c',round=True)
         polythres=10            
    elif region == 'nps2': # As nps, but rectagle like NSIDC # nsidc.org/data/polar_stereo/ps_grids.html
         mapparams = dict(projection='stere',llcrnrlat=33.92,urcrnrlat=31.37,
                          llcrnrlon=279.26,urcrnrlon=102.34,lat_0=90,lon_0=135,
                          resolution='c')
    elif region == 'nps3': #As nps2, but zoomed in
         lat1=37; lat2=36.5;lon1=280;lon2=105;lat_0=90;lon_0=138                       
         mapparams = dict(projection='stere',llcrnrlat=lat1,urcrnrlat=lat2,
                         llcrnrlon=lon1,urcrnrlon=lon2,lat_0=lat_0,lon_0=lon_0,
                          resolution='c')
    elif region == 'baf': #As nps2, but zoomed in
         lat1=50; lat2=80;lon1=-100; lon2=-40
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')

    ##########Subregions
    elif region == 'nh': # NH, standard (cylindrical) projections          
         lat1=0; lat2=90;lon1=0; lon2=360
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif region == 'tpo_na': #Tropical Pacific Ocean + North America (cylindrical)          
         lat1=-20; lat2=80;lon1=100; lon2=310
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif region == 'atl': #Tropical Pacific Ocean + North America (cylindrical)          
         lat1=-30; lat2=80;lon1=-100; lon2=10
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif region == 'na': #North Atlantic          
         lat1=0; lat2=80;lon1=-70; lon2=10
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    else:
        print "Incorrect region. Global: glob, glob_rob"
        print "                  Polar: nps, nps2 nps3 sps"
        print "             Subregions: nh tpo_na"
        
        
        return -1

    mapparams['ax'] = ax
    #mapparams['area_thresh']=100000
    #### make map        
    bm=Basemap(**mapparams)

    #### fill colors        
#    landfillcol=[0.85,0.85,0.85] 
#    oceanfillcol=np.array([176, 237, 245])/255.

    if landmask:
        if oceanmask:
            bm.fillcontinents(landfillcol,lake_color=oceanfillcol)
            bm.drawmapboundary(fill_color=oceanfillcol)
        else:
            bm.fillcontinents(landfillcol)
           
    #### Remove annoying rivers: http://stackoverflow.com/questions/14280312/world-map-without-rivers-with-matplotlib-basemap   
    coasts = bm.drawcoastlines(zorder=1,color='white',linewidth=0) #remove standard coastlines and get handles
    coasts_paths = coasts.get_paths()
    ipolygons = np.arange(polythres)

    for ipoly in range(len(coasts_paths)):
        r = coasts_paths[ipoly]
        # Convert into lon/lat vertices
        polygon_vertices = [(vertex[0],vertex[1]) for (vertex,code) in 
                                    r.iter_segments(simplify=False)]
        px = [polygon_vertices[i][0] for i in xrange(len(polygon_vertices))]
        py = [polygon_vertices[i][1] for i in xrange(len(polygon_vertices))]
        if ipoly in ipolygons:
            bm.plot(px,py,linewidth=coastlinewidth,zorder=3,color='black') # plot only larger lakes
        elif landmask: 
            bm.plot(px,py,linewidth=coastlinewidth,zorder=4,color=landfillcol) # if landmask: fill in with landfillcol
                                                                    # if no land mask: donot plot, otherwise
                                                                    # this will show up on contourf plots
                

    #### drawgrid        
    if drawgrid:
        bm.drawparallels(np.arange(-90.,90.,30.),linewidth=0.2)
        #bm.drawmeridians(np.arange(0.,360.,30.),linewidth=0.2)
    return bm

############################################################################
def make_vax(ax,
             region='glob',plevtop=10, force_xaxis=False,plog='True'):


    """Make vertical axis for lat,plev plots 

    Parameters:
    -----------
      ax : [axis]
            Subplot axis

      *region* : ['glob' | string] 
                Region for the plots. Current options:
                'glob': [-90,90] 
                'nh': [0,90] 

      *plevtop* : [ *10.* | float ]         
                Top pressure level

      Returns:
      --------
          Nothing
    """      
    if ax is None:
       ax=plt.gca()  
    if region=='NH':
        latmin=0
        latmax=85
    elif region=='glob':
        latmin=-85
        latmax=85

    # X-axis
    ax.set_xlim(latmin,latmax)
    if region=='NH':
        ax.set_xticks([0, 20, 40,60, 80])

    if ax.is_last_row() or force_xaxis: 
      if region=='NH':
        ax.set_xticklabels(('EQ', '20N', '40N', '60N','80N'))

    else:
        ax.set_xticklabels((''))
    # Y-axis
    if plog: ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_ylim(1000,plevtop)
    #if plevtop==10:
    ax.set_yticks([1000,500, 200, 100, 50, 20,10])
    if ax.is_first_col(): 
       ax.set_yticklabels((1000,500, 200, 100, 50,20,10))  
       ax.set_ylabel('Pressure')            
    else: ax.set_yticklabels((''))
    # X+Y-axis
    ax.tick_params(which = 'both', direction = 'out')
    ax.minorticks_off() 
   


############################################################################
#COLORPLOTS#################################################################
############################################################################

############################################################################
def add_cf(ax,x,y,fld,
                 clevs=None,
                 cint=None,cint0=None,coffset=0.,nclevspos=5,
                 cmap='blue0red19',cmapi=None,
                 plot_co=True,
                 latlon=True): 

    """add contourf plot to axis.

    Parameters:
    -----------
      ax : [axis]
            Either the basemap (latlon) or figure axis (latpres)

      x : [array] 
            A 1D array with x-dimensions (lon or lat) 

      y : [array] 
            A 1D array with y-dimensions (lat or plev)        

      fld : [array] 
            A 2D array with the data to be plotted

      *clevs* : [ *None* | array or list]  
                A 1D array or list with clevs

      *cint0* : [ *None* | float] 
                Contour interval value used to create clevs centred around 0

      *cint* : [ *None* | float] 
                Contour interval value used to create clevs that (includes 0) 

      *coffset* : [ *0* | float] 
                  Offset for creating clevs from cint        
                  
      *nclevspos* : [ *5* | int] 
                    Number of positive contour levels, used to Contour interval value 
                    used to create clevs that include 0 
             

      *cmap* : [ *'blue0red19'* | string ]   
                        blue0red19: default
                        bluegrayred19: 
                        blue0red11giss: Mimicing GISS T-trends (used for hiatus studies)
                        bluegrayred9dark:  Used for retreat/advance skill studies
                        bluegrayred11dark: Used for retreat/advance skill studies          

      *cmapi* [ *None* | list ]
               list of color indices
      
      *plot_co*: [ *True | boolean ]
               If True (default) plot contour

      Returns:
      --------
          cf: contourf handle (typically to add colorbar after func call)
          
      Notes:
      --------    
          Contourf is more flexible than pcolormesh in the sense that it can 
          automap colors once you specify clevs. 
          Therefore, it supports nclevspos if cint is specified            
          
    """      

    ##### Determine if latlon, and make grid
    #if np.max(x)>90 and np.max(y)<359: #X is lon (not lat), Y is lat (not pres)  
    if latlon==True:    
      # Add cyclic lon if needed
      if np.mod(x.shape,2) == 0:
        fld,x = mpltk.basemap.addcyclic(fld,x)
      if np.max(y)<90: 
            fld,y =add_poles(fld,y)

    xs,ys=np.meshgrid(x,y)        
    cfparams={'latlon':latlon}


    #############clevs and colors################################### 
    cfparams['extend']='both'

    if cmapi is not None: # clevs have to be given, manual mapping
        if clevs is None:
            print 'Error: when cmapi defined, clevs have to be defined too'
            return -1
        elif cint is not None or cint0 is not None:
            print 'Error: when cmapi defined, cint/cint0 cannot be defined'
            return -1
        else:
            colors=plt.get_cmap(cmap)(cmapi)
            cmap, norm = from_levels_and_colors(clevs, colors ,extend='both')
            cfparams['cmap']=cmap           
            cfparams['norm']=norm       
            cfparams['levels']=clevs
    else: # automatically mapping
        cfparams['cmap']=cmap           
        if clevs is not None:
            cfparams['levels']=clevs
        if cint0 is not None: #create clevs centred around 0
            crange=(nclevspos-0.5)*cint0; 
            cmin=-crange+coffset; cmax=crange+coffset
            cfparams['levels']=np.arange(cmin,cmax+cint0,cint0);
        if cint is not None:#create clevs that includes 0
            crange=nclevspos*cint; 
            cmin=-crange+coffset; cmax=crange+coffset
            cfparams['levels']=np.arange(cmin,cmax+cint,cint);
    #cfparams['zorder']=2   
    #############plot######################    
    cf=ax.contourf(xs,ys,fld,**cfparams)
    if plot_co:
       add_co(ax,x,y,fld,
              clevs=clevs,cint=cint,cint0=cint0,coffset=coffset,
              color='gray',linewidth=0.2,neg_dash=False)             
    return cf


############################################################################
def add_pc(ax,x,y,fld,
               clevs=None,
               cint=None,cint0=None,coffset=0.,
               cmap='blue0red19',cmapi=None,extend='both'): 
    """add pcolormesh to axis.
    
 
    Parameters:
    -----------
      ax : axis
            Either the basemap (latlon) or figure axis (latpres)

      x : array 
            A 1D array with x-dimensions (lon or lat) 

      y : array 
            A 1D array with y-dimensions (lat or plev)        

      fld : array 
            A 2D array with the data to be plotted

      *clevs* : [ *None* | array or list]  
                A 1D array or list with clevs

      *cint0* : [ *None* | float] 
                Contour interval value used to create clevs centred around 0 

      *cint* : [ *None* | float] 
                Contour interval value used to create clevs (includes 0) 

      *coffset* : [ *0* | float] 
                  Offset for creating clevs from cint         

      *cmap* : [ *'blue0red19'* | string ]   
                        blue0red19: default
                        bluegrayred19: 
                        blue0red11giss: Mimicing GISS T-trends (used for hiatus studies)
                        bluegrayred9dark:  Used for retreat/advance skill studies
                        bluegrayred11dark: Used for retreat/advance skill studies          

      *cmapi* [ *None* | list ]
               list of color indices
               
      *extend* [ *both* | neither ]
               Colorbar that extends beyond upper limit
               
               

      Returns:
      --------
          pc: pcolormesh handle (typically to add colorbar after func call)

      Notes:
      --------    
          For pcolormesh, clevs not given explicitely, but through 'norm'
          Therefore, if clevs/cint/cint0 given, we have to manually map
          In case of specified cint/cint0, the range is determined by the ncolor in the cmap
               
      Known issues/bugs:
      --------
          1) For lat,pres plot, have to shift the data
          2) Error if basemap is 'glob'

    """      

    ##### Determine if latlon, and make grid
    if np.max(x)>90 and np.max(y)<359: #X is lon (not lat), Y is lat (not pres)  
        # Add cyclic lon if needed
        if np.mod(x.shape,2) == 0:
            fld,x = mpltk.basemap.addcyclic(fld,x)
        lon_edge,lat_edge=centre_to_edge(x,y)
        xs, ys = np.meshgrid(lon_edge,lat_edge) ###NB: len(lat_edge)=nlat+1
        latlon=True
    else: #lat,plev plot
        latlon=False
        xs,ys=np.meshgrid(x,y)
    fld = np.ma.masked_invalid(fld) #http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values                


   #############colors #########
    if cmapi is not None: # create colors from cmapi
        colors=make_colors(cmap,cmapi)
    else: # chose all colors from cmapi
        colors=make_colors(cmap)
        if cint is not None: # remove the middle color if clevs include 0
            ncolors=np.shape(colors)[0]
            ncolorsmin=int(np.floor(ncolors/2.))  
            colors=make_colors(cmap,range(ncolorsmin)+range(ncolorsmin+1,ncolors))

   #############clevs######################################
 
    ncolors=np.shape(colors)[0]
    if clevs is not None:
        clevs=np.asarray(clevs)
    if cint:
        crange=cint*(ncolors/2.-1); 
        clevs=np.linspace(-crange,crange,num=ncolors-1,endpoint=True)+coffset; 

    if cint0:
        crange=cint0*(ncolors/2.-1); 
        clevs=np.linspace(-crange,crange,num=ncolors-1,endpoint=True)+coffset; 
 
    ###########color mapping#####################################

    if clevs is not None:   
        cmap, norm = from_levels_and_colors(clevs, colors ,extend=extend)
        pcparams=dict(cmap=cmap,norm=norm)           
    else:
        pcparams=dict(cmap=cmap)           
        
    if latlon: pcparams['latlon']=latlon 
    #############plot######################    
    pc=ax.pcolormesh(xs,ys,fld,**pcparams)

    return pc
    
############################################################################

def add_cb(ax,pc,
             units=None,
             x0scale=1.,y0scale=1.,lscale=1.,labelsize=10,
             manticks=None,manlabels=None,
             orientation='vertical',
             spacing='proportional'):

    """Adds a cbar to a plot in its own axis
    
    Parameters:
    ----------
      ax : [axis]
            Subplot axis

      pc :  
          pcolormesh or contourf handle

      *x0scale*: [*1.*, float]
          x0 scaling factor

      *y0scale*: [*1.*, float]
          y0 scaling factor

      *lscale*: [*1.*, float]
          length scaling factor (height for vertical, width for horizontal colorbars)
          
      *labelsize*: [*10.*, float]
          height scaling factor

      *manticks*:[None]
          list of colorbar ticks  

      *manlabels*:[None]                                                                                                                                  
          list of colorbar labels                                                                                                       

      *orientation*:[*'vertical'*,string]                                                                                                                                  
          orientation of colorbar

      *spacing*:[*'proportional'*,string]                                                                                                                                  
          spacing of colorbar


    Returns:
    -----------
        Nothing    
    """

    #make axis 
    box = ax.get_position()
    fig = ax.get_figure()
    if orientation=='vertical':
      cbar_ax=fig.add_axes([box.x0+box.width*1.035*x0scale, box.y0*y0scale, 0.02, box.height*lscale])
    if orientation=='horizontal':
      cbar_ax=fig.add_axes([box.x0+box.width*(0.08+x0scale-1), box.y0+box.height*0.08*y0scale, box.width*lscale, 0.015])
      
      
    #plot cbar 
    cbar=fig.colorbar(pc, cax=cbar_ax,extendfrac='auto',orientation=orientation,spacing=spacing)
    cbar.ax.tick_params(labelsize=labelsize) 
    #units 
    if units is not None: cbar.set_label(units)
    if manticks is not None:
        cbar.set_ticks(manticks) 
        cbar.set_ticklabels(manticks)
        if manlabels is not None:
           cbar.set_ticklabels(manlabels)
############################################################################
#CONTOURPLOTS#################################################################
############################################################################

############################################################################
def add_co(ax,x,y,fld,
                clevs=None,cint=None,cint0=None,coffset=0.,nclevspos=50,
                color='darkgray',linewidth=1.5,
                neg_dash='True',clevs_color=None):
    """add contour plot to axis.
    
    Parameters:
    -----------
      ax : [axis]
            Either the basemap (latlon) or figure axis (latpres)

      x : [array] 
            A 1D array with x-dimensions (lon or lat) 

      y : [array] 
            A 1D array with y-dimensions (lat or plev)        

      fld : [array] 
            A 2D array with the data to be plotted

      *clevs* : [ *None* | array or list]  
                A 1D array or list with clevs

      *cint0* : [ *None* | float] 
                Contour interval value used to create clevs (centred around 0)

      *cint* : [ *None* | float] 
                Contour interval value used to create clevs that (includes 0) 

      *coffset* : [ *0* | float] 
                  Offset for creating clevs from cint         

      *nclevspos* : [ *50* | int] 
                    Number of positive contour levels, used to Contour interval value 
                    used to create clevs that include 0 
 
      *color* : [ *'darkgray'* | string ]   
                Color of contours 
                
      *linewidth*: [ *1.5* | int]
                   Linewidth of the contours 

      *neg_dash*: [ *True* | boolean ]
                   If true negative contours are dashed       

      *clevs_color*: [ *None* | list]
                     If specified, it will add up to 3 colorcontours (green,blue,red)

      Returns:
      --------
          Nothing

    """      
    ##### Determine if latlon, and make grid
    if np.max(x)>90 and np.max(y)<359: #X is lon (not lat), Y is lat (not pres)  
        latlon=True    
        # Add cyclic lon if needed
        if np.mod(x.shape,2) == 0:
            fld,x = mpltk.basemap.addcyclic(fld,x)
        if np.max(y)<90: 
            fld,y =add_poles(fld,y)
   
    else:
        latlon=False        
    xs,ys=np.meshgrid(x,y)
    coparams=dict(latlon=latlon)  

    #############clevs###################################### 
    if clevs is not None: clevs=np.array(clevs) #convert to array if it was a list
    if cint0 is not None: #create clevs centred around 0
        crange=(nclevspos-0.5)*cint0; 
        cmin=-crange; cmax=crange
        clevs=np.arange(cmin,cmax+cint0,cint0)+coffset;
    if cint is not None: #create clevs that includes 0
        crange=nclevspos*cint; 
        cmin=-crange; cmax=crange
        clevs=np.arange(cmin,cmax+cint,cint)+coffset;

    #############contour settingss############################### 
    coparams['colors']=color
    coparams['linewidths']=linewidth
    #############Plot############################### 

    if clevs is None:
        ax.contour(xs,ys,fld,**coparams)
    else:   
        if np.max(clevs)>0:
            coparams['levels']=clevs[clevs>0.]  
            ax.contour(xs,ys,fld,**coparams)
        if np.min(clevs)<0:
            coparams['levels']=clevs[clevs<0.]
            CS_neg=ax.contour(xs,ys,fld,**coparams)        
            
            for c in CS_neg.collections:
                if neg_dash:
                    c.set_dashes([(0, (3.0, 2.0))])
                else:
                    c.set_dashes((None,None))            
        if 0. in clevs:
            coparams['levels']=[0]  
            coparams['linewidths']=linewidth*2              
            ax.contour(xs,ys,fld,**coparams)
            
    #additional contours
    if clevs_color is not None:  
        if len(clevs_color)>=1:
           ax.contour(xs,ys,fld,levels=[clevs_color[0]],colors='g',linewidths=1,latlon=latlon)
        if len(clevs_color)>=2:
           ax.contour(xs,ys,fld,levels=[clevs_color[1]],colors='b',linewidths=1,latlon=latlon)
        if len(clevs_color)>=3:
           ax.contour(xs,ys,fld,levels=[clevs_color[2]],colors='r',linewidths=1,latlon=latlon)

############################################################################
#SCATTERPLOTS#################################################################
############################################################################

############################################################################
def add_sc(ax,x,y,fld):
    """add scatter plot to axis.
    
    Parameters:
    -----------
      ax : [axis]
            Either the basemap (latlon) or figure axis (latpres)

      x : [array] 
            A 1D array with x-dimensions (lon or lat) 

      y : [array] 
            A 1D array with y-dimensions (lat or plev)        

      fld : [array] 
            A 2D array with the data to be plotted


"""
    ##### Determine if latlon, and make grid
    if np.max(x)>90 and np.max(y)<359: #X is lon (not lat), Y is lat (not pres)  
        latlon=True    
        # Add cyclic lon if needed
        if np.mod(x.shape,2) == 0:
            fld,x = mpltk.basemap.addcyclic(fld,x)
        if np.max(y)<90: 
            fld,y =add_poles(fld,y)
    else:
        latlon=False
    xs,ys=np.meshgrid(x,y)
    print latlon
    ax.scatter(xs,ys,s=fld)      

############################################################################
def add_sc2d(ax,x,y):
    """2d scatterplot.
    
    Parameters:
    -----------
      ax : [axis]
            Either the basemap (latlon) or figure axis (latpres)

      x : [array] 
            A 1D array with x-data 

      y : [array] 
            A 1D array with y-data         
     """

    ax.scatter(x,y)
    



############################################################################
#UTILS#######################################################################
############################################################################

def add_poles(fld,lat):

    nlat=np.shape(fld)[0];
    nlon=np.shape(fld)[1];

    lat_new=np.zeros(nlat+2)
    fld_new=np.ma.zeros((nlat+2,nlon))

    if(lat[1]>lat[0]):
        lat_new[0]=-90
        lat_new[1:-1]=lat
        lat_new[-1]=90

    fld_new[0,:]=np.mean(fld[0,:])
    fld_new[1:-1,:]=fld
    
    for ilon in range(nlon): fld_new[-1,ilon]=np.ma.masked

    return fld_new,lat_new
        

##save_fig##########################################################################
def mysavefig(fig,fnamesave):
    fig.savefig(fnamesave,dpi=400,bbox_inches='tight')


##make_colors##########################################################################

def make_colors(cmap_name,cmapi=None):
    """Make array of colors rgb values 

    Returns:
    -----------
      Colors: array
        Array of rgb values of colors

    Notes:
    -----------
        Only to be used in rms_plots.py
    """
    cmap=plt.get_cmap(cmap_name)
    if cmapi is None:
        cmapi=np.arange(cmap.N)
    colors=cmap(np.asarray(cmapi))
    #return 
    return colors



##qplot##########################################################################

def qplot(x,y,fld,
            fnamesave=None):
                
    """Quick plot

            
    Parameters:
    -----------
      x : [array] 
            A 1D array with x-dimensions (lon or lat) 

      y : [array] 
            A 1D array with y-dimensions (lat or plev)        

      fld : [array] 
            A 2D array with the data to be plotted

      *fnamesave* [ *None* | string]
            If given save the plot into the file fnamesave

    Returns:
    -----------
        Figure handle    
    """
        
    #creat figure and axis
    fig=plt.figure()    
    ax=plt.gca()
    #determine latlon
    ##### Determine if latlon, and make grid
    if np.max(x)>90 and np.max(y)<359: #X is lon (not lat), Y is lat (not pres)  
        latlon=True    
    else:    
        latlon=False    

    #create map and plot
    if latlon:
        bm=make_bm(ax=ax)
        pc=add_cf(bm,x,y,fld)
    else:
        make_vax(ax=ax)
        pc=add_cf(ax,x,y,fld)
    #cbar
    add_cb(ax,pc)
    #save 
    if fnamesave is not None:
        mysavefig(fig,fnamesave)
    #return 
    return fig



############################################################################

def add_title(ax,title,
              panellab=None,fontsizelab=14):
    """add title to axis.
    
    Parameters:
    -----------
      ax : [axis]
            Subplot axis

      title: [string]
             Title

      *panellab*: [ 'None' | string ]
                  Panel label (typically a,b,c,d..)

      *fontsizelab*: [ '14', float ]
                Fontsize of label 

      Returns:
      --------
          Nothing     
"""
    ax.set_title(title)
    if panellab is not None:
        ax.annotate(panellab,xy=(0.01,1.02),
                          xycoords='axes fraction',fontsize=fontsizelab,fontweight='bold')
      

############################################################################
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

############################################################################
    
############################################################################
#ABOLETE####################################################################
############################################################################

#basemap (for latlon plots)##################################################
def make_basemap(maptype='glob',drawgrid=False,axis=None,
                 landfillcol=[0.85,0.85,0.85],
                 oceanfillcol=np.array([176, 237, 245])/255.,
                 coastlinewidth=0.1,latps0=50,lonps0=270):

    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: make_hori_plot_axis
    """      

    if maptype == 'rob':
       mapparams = dict(projection='robin',lon_0=0,resolution='l')
    elif maptype == 'nps': #North Pole stereographic 
         mapparams = dict(projection='npstere',boundinglat=latps0,lon_0=lonps0,
                     resolution='c',round=True)
    elif maptype == 'nsnh': # NSIDC-Northern Hemispher
                             # https://nsidc.org/data/polar_stereo/ps_grids.html
         mapparams = dict(projection='stere',llcrnrlat=33.92,urcrnrlat=31.37,
                          llcrnrlon=279.26,urcrnrlon=102.34,lat_0=90,lon_0=135,
                          resolution='c')
    elif maptype == 'nsnhz': #zoomed in version of nsidc map
         lat1=37; lat2=36.5;lon1=280;lon2=105;lat_0=90;lon_0=138                       
         mapparams = dict(projection='stere',llcrnrlat=lat1,urcrnrlat=lat2,
                         llcrnrlon=lon1,urcrnrlon=lon2,lat_0=lat_0,lon_0=lon_0,
                          resolution='c')
    elif maptype == 'glob': #latlon          
         lat1=-90; lat2=90;lon1=-180; lon2=180
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif maptype == 'tpo_na': #latlon          
         lat1=-20; lat2=80;lon1=100; lon2=310
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    elif maptype == 'nh': #latlon          
         lat1=0; lat2=90;lon1=0; lon2=360
         mapparams = dict(projection='cyl',llcrnrlat=lat1,urcrnrlat=lat2,\
            llcrnrlon=lon1,urcrnrlon=lon2,resolution='c')
    else:
        print "Incorrect maptype. glob,nps, nsnh, or nsnhz"
        return -1

    if axis != None: # if an axis is given, add to dict for basemap
        mapparams['ax'] = axis
                     
    map=Basemap(**mapparams)
    map.drawcoastlines(linewidth=coastlinewidth,color="0.1")
    if landfillcol != '':
        if oceanfillcol !='':
           map.fillcontinents(landfillcol,lake_color=oceanfillcol)
           map.drawmapboundary(fill_color=oceanfillcol)
        else:
           map.fillcontinents(landfillcol)

    #### drawgrid        
    if drawgrid:
        map.drawparallels(np.arange(-90.,90.,30.),linewidth=0.2)
        map.drawmeridians(np.arange(0.,360.,30.),linewidth=0.2)
    return map



   
############################################################################
def add_plot_latlon(lat,lon,fld,map,axis=None,cint='',clevs='',colors='',
                shadetype='contourf',supp_cb='False',extend='both'): 

    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: add_colorplot_latlon
    """      
    ##### Build grid, add cyclic lon if needed, shift lons+lats in case of pcolor plot
    if np.mod(lon.shape,2) == 0:
       fld,lon = mpltk.basemap.addcyclic(fld,lon)
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
        pc=map.pcolormesh(lons,lats,fld,**pcparams)
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

############################################################################
def add_color_hori(ax,lon,lat,fld,cint=None,clevs=None,colors=None,
                shadetype='contourf',supp_cb='False',extend='both'): 

    """
    NOTE: Old name: add_colorplot_latlon
    """      
    ##### Build grid, add cyclic lon if needed, shift lons+lats in case of pcolor plot
    if np.mod(lon.shape,2) == 0:
       fld,lon = mpltk.basemap.addcyclic(fld,lon)
    if shadetype == 'pcolor':
        lon_edge,lat_edge=centre_to_edge(lon,lat)
        lons, lats = np.meshgrid(lon_edge,lat_edge) ###NB: len(lat_edge)=nlat+1
        fld = np.ma.masked_invalid(fld) #http://stackoverflow.com/questions/7778343/pcolormesh-with-missing-values                
    else:
        lons, lats = np.meshgrid(lon,lat)  
    #############Colors
    if colors is None:
        colors=plt.get_cmap('bluegrayred19')(np.arange(19))

    #############clevs###################################### 
    # options:
    # 1) give clevs --> clevs=clevs
    # 2) give cint --> clevs are constructed centred around 0
    # 3) no clevs --> 19 autolevels

    ncolors=colors.shape[0]; #nr of colors in cbar
    if cint:  #construct from cint (if given)
       crange=cint*(ncolors-2)/2; cmin=-crange; cmax=crange
       clevs=np.linspace(cmin,cmax,ncolors-1);
    #############CMAP######################################
    if clevs:        
       cmaploc, norm = from_levels_and_colors(clevs, colors ,extend=extend)
       pcparams=dict(latlon=True,cmap=cmaploc,norm=norm)
    else:
       cmaploc=col.ListedColormap(colors)         
       pcparams=dict(latlon=True,cmap=cmaploc)
    #############plot######################    
    if shadetype == 'pcolor':
       pc=ax.pcolormesh(lons,lats,fld,**pcparams)
    elif shadetype == 'contourf':
         if clevs:
             pcparams['levels']=clevs
             pcparams['extend']='both'
             pc=ax.contourf(lons,lats,fld,**pcparams)
         else:
             pc=ax.contourf(lons,lats,fld,19,**pcparams)    
    else:
        print "Incorrect shadetype. Choose pcolor or contourf"
        return -1
    ############cbar######################
    if supp_cb==False:
        if clevs:
            ax.colorbar(pc,location='right',pad="5%",ticks=clevs)
        else:
            ax.colorbar(pc,location='right',pad="5%")
    return pc

############################################################################
def add_color_vert(ax,lon,lat,fld,cint=None,clevs=None,colors=None,
                shadetype='contourf',vmin=None,vmax=None,supp_cb='False',extend='both'): 

    """
    NOTE: Old name: add_colorplot_latlon
    """      
    ##### Build grid, add cyclic lon if needed, shift lons+lats in case of pcolor plot
    lons, lats = np.meshgrid(lon,lat)  
    #############Colors
    if colors is None:
        colors=plt.get_cmap('bluegrayred19')(np.arange(19))

    #############clevs###################################### 
    # options:
    # 1) give clevs --> clevs=clevs
    # 2) give cint --> clevs are constructed centred around 0
    # 3) no clevs --> 19 autolevels

    ncolors=colors.shape[0]; #nr of colors in cbar
    if cint:  #construct from cint (if given)
       crange=cint*(ncolors-2)/2; cmin=-crange; cmax=crange
       clevs=np.linspace(cmin,cmax,ncolors-1);
    #############CMAP######################################
    if clevs:        
       cmaploc, norm = from_levels_and_colors(clevs, colors ,extend=extend)
       pcparams=dict(cmap=cmaploc,norm=norm)
    else:
       cmaploc=col.ListedColormap(colors)         
       pcparams=dict(cmap=cmaploc)
    #############plot######################    
    if shadetype == 'pcolor':
        print pcparams
        pc=ax.pcolormesh(lons,lats,fld,**pcparams)       
    elif shadetype == 'contourf':
         if clevs != '':
             pcparams['levels']=clevs
             pcparams['extend']='both'
             pc=ax.contourf(lons,lats,fld,**pcparams)
         else:
             pc=ax.contourf(lons,lats,fld,19,**pcparams)    
    else:
        print "Incorrect shadetype. Choose pcolor or contourf"
        return -1
    ############cbar######################
    if supp_cb==False:
        if clevs:
            ax.colorbar(pc,location='right',pad="5%",ticks=clevs)
        else:
            ax.colorbar(pc,location='right',pad="5%")
    return pc



############################################################################
def make_vert_plot_axis(ax,title,region='NH',plevtop=10):
    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: make_vax
    """      
    if region=='NH':
        latmin=0
        latmax=85
    # X-axis
    ax.set_xlim(latmin,latmax)
    if region=='NH':
        ax.set_xticks([0, 20, 40,60, 80])
    if ax.is_last_row(): 
        if region=='NH':        
            ax.set_xticklabels(('EQ', '20N', '40N', '60N','80N'),fontsize=6,fontweight='bold')
    else:
        ax.set_xticklabels((''))
    # Y-axis
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_ylim(1000,plevtop)
    if plevtop==10:
        ax.set_yticks([1000,500, 200, 100, 50, 20,10])
        if ax.is_first_col(): 
            ax.set_yticklabels((1000,500, 200, 100, 50,20,10),fontsize=7,fontweight='bold')  
            ax.set_ylabel('Pressure',fontsize=9)            
        else:
            ax.set_yticklabels((''))
    # X+Y-axis
    ax.tick_params(which = 'both', direction = 'out')
    ax.minorticks_off() 
    # Title
    ax.set_title(title,fontsize=8)    

############################################################################

def add_vert_plot_contour(ax,fld,plev,lat,clevs='',cint0='',cint='',clevs_color='',label_clabs=False):
    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: add_contour
    """      

    # prep
    lats,levs=np.meshgrid(lat,plev)
    #clevs
    # options:
    # 1) give clevs --> clevs=clevs
    # 2a) give cint0 --> clevs are constructed centred around 0
    # 2b) give cint  --> clevs that includes 0
    # 3) no clevs --> 19 autolevels
    nclevs=50

    if cint0 !='':  #construct from cint0 (if given)
       crange=cint0*(nclevs-1)/2.; cmin=-crange; cmax=crange
       clevs=np.arange(cmin,cmax+cint0,cint0);
    if cint !='':  #construct from cint (if given)
       crange=cint*nclevs/2.; cmin=-crange; cmax=crange
       clevs=np.arange(cmin,cmax+cint,cint);

    clevs_neg=clevs[clevs<0.]
    clevs_pos=clevs[clevs>0.]

    #plot
    ax.contour(lats,levs,fld,levels=clevs_pos,colors='k',linewidths=0.5)
    CS_neg=ax.contour(lats,levs,fld,levels=clevs_neg,colors='k',linewidths=0.5)        
    for c in CS_neg.collections:
        c.set_dashes([(0, (3.0, 2.0))])
    #additional contours
    if 0. in clevs:
       ax.contour(lats,levs,fld,levels=[0],colors='k',linewidths=1)
    if clevs_color !='':  
        if len(clevs_color)>=1:
           ax.contour(lats,levs,fld,levels=[clevs_color[0]],colors='g',linewidths=1)
        if len(clevs_color)>=2:
           ax.contour(lats,levs,fld,levels=[clevs_color[1]],colors='b',linewidths=1)
        if len(clevs_color)>=3:
           ax.contour(lats,levs,fld,levels=[clevs_color[2]],colors='r',linewidths=1)
      
############################################################################

def add_vert_plot_contourf(ax,fld,plev,lat,colors='',clevs='',cint='',cint0=''):
    """
    NOTE: This routine is now absolete, but was frequently used prior to 20161123
          New name: add_color_vert
    """      


    # prep
    lats,levs=np.meshgrid(lat,plev)
    #colors 
    if colors=='':
       full_cmap = plt.get_cmap('bluered19')
       colors=full_cmap(np.arange(19))
       nclevs=colors.shape[0]-1; #nr of colors in cbar - 1
    if cint0 !='':  #construct from cint0 (if given) [centred around 0]
       crange=cint0*(nclevs-1.)/2.; cmin=-crange; cmax=crange
       clevs=np.arange(cmin,cmax+cint0,cint0);
    if cint !='':  #construct from cint (if given)# [not centred around 0 -->remove the middle color]
       crange=cint*nclevs/2.; cmin=-crange; cmax=crange
       clevs=np.arange(cmin,cmax+cint,cint);
    ax.contourf(lats,levs,fld,levels=clevs,colors=colors) 

############################################################################
def plot_latlon(lat,lon,fld,axis=None,fnamesave='',
                maptype='glob',drawgrid=False,landfillcol='',oceanfillcol='',coastlinewidth=0.1,
                cint='',clevs='',colors='',shadetype='contourf',supp_cb=False,extend='both',
                title='',panellab=''):
    """
    plot_latlon(lat,lon,fld,axis=None,cint='',clevs='',
                colors='',shadetype='contourf',maptype='glob',
                supp_cb=False,drawgrid=False,title='',panellab='',fnamesave='')
    Inputs: lat: 1D lat array
            lon: 1D lon array
            fld: 2D matrix of data [lat x lon]
            axis: handle to the axes instance to give to basemap. default: None
            colors: colors
            clevs: contour levels (can also be controlled by cint)
            cint: contourinterval for clevs centred around 0            
            shadetype: Shading method (contourf or pcolor). default: pcolor
            maptype: Map type (glob or nhpc). default: glob 
            supp_cb: Suppress colorbar. default: False
            drawgrid: if True, draw parallels and meridians. default: False
            title: title
            panellab: label of panel (a,b,c, etc..)              
            fnamesave: If defined, print output to fnamesave.png and fnamesave.pdf
    Returns: basemap handle (to add to bm after function call),            
             pcolormesh handle (typically to add colorbar after func call)
   Examples: Testplot: plot_latlon(lat,lon,fld,fnamesave='test')         
             
             
    """
    ########################################################
    #############PREPS######################################
    ########################################################

    ## open new figure if fnamesave is defined
    if fnamesave !='': fig=plt.figure()
    bm=make_basemap(axis=axis,maptype=maptype,drawgrid=drawgrid,
                   landfillcol=landfillcol,oceanfillcol=oceanfillcol,coastlinewidth=coastlinewidth)
    pc=add_plot_latlon(lat,lon,fld,bm,axis=None,
            shadetype=shadetype,clevs=clevs,colors=colors,supp_cb=supp_cb,extend=extend)
    add_title(title=title,panellab=panellab)                  
     ############save figure######################
    if fnamesave !='': fig.savefig(fnamesave + '.png',dpi=400) 
    return map,pc


##Colormaps - old##########################################################################


def register_rms_cmaps_old():
    """create my personal colormaps with discrete colors and register them.
       default is to register all of them. can also specify which one.
       (@@ input arg cmap not implemented yet 2/27/14)
    """
    #print 'registering cmaps'
###############
    # bluegrayred19
    # blueish at top, gray in middle, reddish at bottom

    colors = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [200, 200, 200],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)/255
    
    thecmap = col.ListedColormap(colors,'bluegrayred19')
    cm.register_cmap(cmap=thecmap)

    # ============================================
    # bluegrayred19_r (the above, flipped)
    #
    bluered19_r = np.flipud(colors)
    thecmap = col.ListedColormap(bluered19_r,'bluegrayred19_r')
    cm.register_cmap(cmap=thecmap)
###############
    # bluered19
    # blueish at top, white in middle, reddish at bottom

    colors = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [255, 255, 255],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)/255
    
    thecmap = col.ListedColormap(colors,'bluered19')
    cm.register_cmap(cmap=thecmap)

    # ============================================
    # blue2red19_r (the above, flipped)
    #
    bluered19_r = np.flipud(colors)
    thecmap = col.ListedColormap(bluered19_r,'bluered19_r')
    cm.register_cmap(cmap=thecmap)

##################
    # bluegrayred9
    # Adapted 11-class RdBu from colorbrewer2.org:
    # Recipe:
    # 1) Pick the 4 darkest and 4 lightest colors from 11-class Rdbu
    # 2) Replace the 3 middle ones with a gray shading (colorblind friendly)
    colors = np.array([ [5,48,97], \
                       [33,102,172], \
                       [67,147,195],\
                       [146,197,222],\
                       [130,130,130],\
                       [244,165,130],\
                       [214,96,77],\
                       [178,24,43],\
                       [103,0,31]],\
                    dtype=float)/255.

    thecmap = col.ListedColormap(colors,'bluegrayred9')
    cm.register_cmap(cmap=thecmap)
    # ============================================
    # bluegrayred9_r (the above, flipped)
    #
    bluegrayred9_r = np.flipud(colors)
    thecmap = col.ListedColormap(bluegrayred9_r,'bluegrayred9_r')
    cm.register_cmap(cmap=thecmap)

##################
    # bluegrayred11
    # Adapted 11-class RdBu from colorbrewer2.org:
    # Recipe:
    # 1) 11-class Rdbu
    # 2) Replace the white color with gray (colorblind friendly)
    colors = np.array([ [5,48,97], \
                       [33,102,172], \
                       [67,147,195],\
                       [146,197,222],\
                       [209,229,240],\
                       [130,130,130],\
                       [253,219,199],\
                       [244,165,130],\
                       [214,96,77],\
                       [178,24,43],\
                       [103,0,31]],\
                    dtype=float)/255.

    thecmap = col.ListedColormap(colors,'bluegrayred11')
    cm.register_cmap(cmap=thecmap)
    # ============================================
    # bluegrayred9_r (the above, flipped)
    #
    bluegrayred11_r = np.flipud(colors)
    thecmap = col.ListedColormap(bluegrayred11_r,'bluegrayred11_r')
    cm.register_cmap(cmap=thecmap)
#
register_rms_cmaps()

