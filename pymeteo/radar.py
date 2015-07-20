#!/usr/bin/env python

import numpy as np
import matplotlib
from pylab import pcolor, colorbar
import os

cmap_radar_level_step = 0.5
"""Used to define the resolution of the color map."""

def cmap_radar_levels_full():
   return list(np.arange(-25,75,cmap_radar_level_step))

def cmap_radar_levels():
   return list(np.arange(0,75,cmap_radar_level_step))

def cmap_radar_colors():
   return ('#ffffff', # -25  0.0   1.0    1.0    1.0     white
           '#00ecec', #   5  0.3   0.0    0.923  0.923   cyan
           '#01A0f6', #  10  0.35  0.004  0.627  0.965   light blue
           '#0000f6', #  15  0.4   0.0    0.0    0.965   blue
           '#00ff00', #  20  0.45  0.0    1.0    0.0     green
           '#00c800', #  25  0.5   0.0    0.784  0.0     green yellow
           '#009000', #  30  0.55  0.0    0.6    0.0     mid green
           '#ffff00', #  35  0.6   1.0    1.0    0.0     yellow
           '#e7c000', #  40  0.65  0.905  0.753  0.0     mid yellow
           '#ff9000', #  45  0.7   1.0    0.656  0.0     dark yellow
           '#ff0000', #  50  0.75  1.0    0.0    0.0     red
           '#d60000', #  55  0.8   0.839  0.0    0.0     mid red
           '#c00000', #  60  0.85  0.753  0.0    0.0     dark red
           '#ff00ff', #  65  0.9   1.0    0.0    1.0     purple
           '#9955c9') #  70  0.95  0.6    0.333  0.788   purpleish 
#           #ffffff   #  75  1.0   1.0    1.0    1.0     white

def cmap_contour_levels():
   return [30, 40, 50, 60, 65]

def cmap_radar():
   cdict = { 'red' : [(0.0,   1.0,   1.0),
                      (0.067, 1.0,   0.0),
                      (0.133, 0.004, 0.004),
                      (0.2,   0.0,   0.0),
                      (0.467, 0.0,   1.0),  
                      (0.533, 0.905, 0.905),
                      (0.6,   1.0,   1.0),
                      (0.667, 1.0,   1.0),
                      (0.733, 0.839, 0.839),
                      (0.8,   0.753, 0.753),
                      (0.867, 0.588, 1.0),
                      (0.933, 0.6,   0.6),
                      (1.0,   1.0,   1.0)], 

           'green':  [(0.0,   1.0,   1.0),
                      (0.067, 0.923, 0.923),
                      (0.133, 0.627, 0.627),
                      (0.2,   0.0,   0.0),
                      (0.267, 1.0,   1.0),
                      (0.333, 0.784, 0.784),
                      (0.4,   0.6,   0.6),
                      (0.467, 0.55,  1.0),
                      (0.533, 0.753, 0.753),
                      (0.6,   0.656, 0.656),
                      (0.667, 0.0,   0.0),
                      (0.933, 0.0,   0.333),
                      (1.0,   0.333, 1.0)],

           'blue' :  [(0.0,   1.0,   1.0),
                      (0.067, 0.923, 0.923), 
                      (0.133, 0.965, 0.965), 
                      (0.2,   0.965, 0.965), 
                      (0.267, 0.0,   0.0), 
                      (0.867, 0.0,   1.0), 
                      (0.933, 0.788, 0.788), 
                      (1.0,   1.0,   1.0)] }

   #return matplotlib.colors.LinearSegmentedColormap('my_radar', cdict)#, lut=512)
   matplotlib.cm.register_cmap(name='pymeteo_radar', data=cdict, lut=512)

def cmap_radar_full():
   cdict = { 'red' : [(0.0, 1.0, 1.0),
                      (0.3, 1.0, 0.0),
                      (0.35, 0.004, 0.004),
                      (0.4, 0.0, 0.0),
                      (0.60, 0.0, 1.0),  
                      (0.65, 0.905, 0.905),
                      (0.7, 1.0, 1.0),
                      (0.75, 1.0, 1.0),
                      (0.8, 0.839, 0.839),
                      (0.85, 0.753, 0.753),
                      (0.9, 0.588, 1.0),
                      (0.95, 0.6, 0.6),
                      (1.0, 1.0, 1.0)], 

           'green':  [(0.0, 1.0, 1.0),
                      (0.3, 0.923, 0.923),
                      (0.35, 0.627, 0.627),
                      (0.4, 0.0, 0.0),
                      (0.45, 1.0, 1.0),
                      (0.5, 0.784, 0.784),
                      (0.55, 0.6, 0.6),
                      (0.6, 0.55, 1.0),
                      (0.65, 0.753, 0.753),
                      (0.7, 0.656, 0.656),
                      (0.75, 0.0, 0.0),
                      (0.95, 0.0, 0.333),
                      (1.0, 0.333, 1.0)],

           'blue' :  [(0.0 , 1.0, 1.0),
                      (0.3,  0.923, 0.923), 
                      (0.35, 0.965, 0.965), 
                      (0.4,  0.965, 0.965), 
                      (0.45, 0.0, 0.0), 
                      (0.9,  0.0, 1.0), 
                      (0.95, 0.788, 0.788), 
                      (1.0 , 1.0, 1.0)] }

   #return matplotlib.colors.LinearSegmentedColormap('my_radar', cdict)#, lut=512)
   matplotlib.cm.register_cmap(name='pymeteo_radar_full', data=cdict, lut=512)


#setup colormaps
cmap_radar()
cmap_radar_full()
