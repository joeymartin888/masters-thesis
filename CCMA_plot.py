#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:25:03 2020

@author: Dr. Michael Sigmond, Canadian Centre for Climate Modelling and Analysis
"""


import matplotlib.colors as col


import matplotlib.cm as cm

import numpy as np


def register_cccmacms(cmap='all'):
    
    
    """create my personal colormaps with discrete colors and register them.
    
    
    default is to register all of them. can also specify which one.
    
    
    (@@ input arg cmap not implemented yet 2/27/14)
    
    
    """
    
    
    #print 'registering cmaps'
    
    
    
    
    
    
    # define individual colors as RGB triples
    
    
    # from colorwheel.m
    
    
    # =============================================
    
    
    # kem_w20 (20) OR blue2red_w20
    
    
    # blueish at top, white in middle, reddish at bottom
    
    
    
    cpool = np.array([ [153,255,255], \
    
    
    [204,255,229], \
    
    
    [240,255,240],\
    
    
    [204,255,153],\
    
    
    [178,255,102],\
    
    
    [216,255,76],\
    
    
    [255,255,51],\
    
    
    [255,220,51],\
    
    
    [255,187,51],\
    
    
    [255,153,51],\
    
    
    [255,0,0],\
    
    
    [204,0,0],\
    
    
    [153,0,0]], \
    
    
    dtype=float)
    
    
    
    acccbar = (cpool/255.)
    
    
    thecmap = col.ListedColormap(acccbar,'acccbar')
    
    
    cm.register_cmap(cmap=thecmap)

    return

register_cccmacms()



