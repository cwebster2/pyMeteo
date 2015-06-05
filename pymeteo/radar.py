#!/usr/bin/env python

import numpy as np
import matplotlib
from pylab import pcolor, colorbar
import os

#translate cmap file to a cmap

def cmap_radar_levels():
   return [-25, 5,10,15,20,25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
   #return [-35, 20, 25, 30, 25, 40, 45, 50, 55, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75]

def cmap_radar_colors():
   return ('#ffffff', 
           '#00ecec',
           '#01A0f6',
           '#0000f6',
           '#00ff00',
           '#00c800',
           '#009000',
           '#ffff00',
           '#e7c000',
           '#ff9000',
           '#ff0000',
           '#d60000',
           '#c00000',
           '#ff00ff',
           '#9955c9')

def cmap_contour_levels():
   return [30, 40, 50, 60, 65]

def cmap_radar():
   cdict = { 'red' : [(0.0 , 1.0  , 1.0),
                      (0.3 , 1.0  , 0.0), 
                      (0.35, 0.004, 0.004), 
                      (0.4 , 0.0, 0.0), 
                      (0.45, 0.0, 0.0), 
                      (0.59 , 0.0, 1.0), 
                      (0.65, 1.0, 0.906), 
                      (0.7 , 0.906, 1.0), 
                      (0.74, 0.1, 1.0), 
                      (0.8 , 1.0, 0.839), 
                      (0.85, 0.839, 0.752), 
                      (0.9 , 0.752, 1.0), 
                      (0.95, 1.0, 0.6), 
                      (1.0, 1.0, 1.0)], 

           'green' : [(0.0 , 1.0, 1.0),
                      (0.3 , 1.0, 0.925), 
                      (0.35, 0.925, 0.627), 
                      (0.4 , 0.627, 0.0), 
                      (0.45, 0.0, 1.0), 
                      #(0.5 , 1.0, 0.784), 
                      #(0.55, 0.784, 0.564), 
                      (0.59 , 0.564, 1.0), 
                      (0.6 , 1.0, 1.0), 
                      (0.65, 1.0, 0.753), 
                      (0.7 , 0.753, 0.564), 
                      (0.74, 0.0, 0.0), 
                      (0.75, 0.0, 0.0), 
                      (0.95, 0.0, 0.333), 
                      (1.0 , 0.333, 1.0)], 

	    'blue' : [(0.0 , 1.0, 1.0),
                      (0.3 , 1.0, 0.925), 
                      (0.35, 0.925, 0.964), 
                      (0.45, 0.964, 0.0), 
                      (0.9 , 0.0, 1.0), 
                      (0.95, 1.0, 0.788), 
                      (1.0 , 0.788, 1.0)] }
   print(cdict)

   #return cdict
   return matplotlib.colors.LinearSegmentedColormap('my_radar', cdict, lut=512)
   #matplotlib.cm.register_cmap(name='my_radar', data=cdict, lut=512)


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

