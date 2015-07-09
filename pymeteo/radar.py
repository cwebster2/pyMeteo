#!/usr/bin/env python

import numpy as np
import matplotlib
from pylab import pcolor, colorbar
import os

#translate cmap file to a cmap

def cmap_radar_levels():
   return [-25, 5,10,15,20,25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75]

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

def cmap_radar_test():
   colors =[(0.0, '#ffffff'),
            (0.3, '#00ecec'),
            (0.35,'#01A0f6'),
            (0.4, '#0000f6'),
            (0.45,'#00ff00'),
            (0.5, '#00c800'),
            (0.55,'#009000'),
            (0.6, '#ffff00'),
            (0.65,'#e7c000'),
            (0.7, '#ff9000'),
            (0.75,'#ff0000'),
            (0.8, '#d60000'),
            (0.85,'#c00000'),
            (0.9, '#ff00ff'),
            (0.95,'#9955c9'),
            (1.0, '#ffffff')]
   return matplotlib.colors.LinearSegmentedColormap.from_list('radar', colors)

def cmap_radar():
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
                      (0.9, 1.0, 1.0),
                      (0.95, 0.6, 0.6),
                      (1.0, 1.0, 1.0)], 

           'green':  [(0.0, 1.0, 1.0),
                      (0.3, 0.923, 0.923),
                      (0.35, 0.627, 0.627),
                      (0.4, 0.0, 0.0),
                      (0.45, 1.0, 1.0),
                      (0.5, 0.784, 0.784),
                      (0.55, 0.6, 0.6),
                      (0.6, 1.0, 1.0),
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

   return matplotlib.colors.LinearSegmentedColormap('my_radar', cdict)#, lut=512)
   #matplotlib.cm.register_cmap(name='pymeteo_radar', data=cdict, lut=512)


def cmap_radar_old():
   print("loading data")

   v,r1,g1,b1,r2,g2,b2,dbz = np.loadtxt(os.path.join(os.path.dirname(__file__),"cmap_radar.txt"), unpack=True)

   n = len(v)
   r1= r1/ 255.
   g1= g1/ 255.
   b1= b1/ 255.
   r2= r2/ 255.
   g2= g2/ 255.
   b2= b2/ 255.
   red = []
   blue = []
   green = []
   for i in range(n):
      red.append((v[i],r1[i],r2[i]))
      green.append((v[i],g1[i],g2[i]))
      blue.append((v[i],b1[i],b2[i]))

   red = tuple(red)
   green = tuple(green)
   blue = tuple(blue)

   cdict = { 'red' : red, 'green' : green, 'blue' : blue }
   print(cdict)

   return matplotlib.colors.LinearSegmentedColormap('radar', cdict, 512)


def main():
   print('in main')
   my_cmap = cmap_radar()

if __name__ == "__main__":
   print('Started standalon, running main')
   main()

