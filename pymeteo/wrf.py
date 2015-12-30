# wrf.py
import numpy as np

######
# The following code is adapted from NCL which is copyright:
#
# Copyright © 2007-2015 University Corporation for Atmospheric Research (UCAR)
#
# and subject to the following license
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# Neither the names of NCAR's Computational and Information Systems Laboratory,
# the University Corporation for Atmospheric Research, nor the names of its
# contributors may be used to endorse or promote products derived from this
# Software without specific prior written permission.  Redistributions of source
# code must retain the above copyright notices, this list of conditions, and the
# disclaimer below.  Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions, and the disclaimer below in the
# documentation and/or other materials provided with the distribution.  THIS
#
# SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
#
# https://www.ncl.ucar.edu/Download/NCL_source_license.shtml

def ll_to_ij(map_proj, truelat1, truelat2, stand_lon, dx, dy, ref_lat, ref_lon, lat, lon):
    
    if (map_proj == 6):
        #TODO: implement map_proj 6
        #pole_lat = f.getncattr('POLE_LAT')
        #pole_lon = f.getncattr('POLE_LON')
        latinc = ( dy * 360. ) / 2.0 / np.pi / 6.37e6
        loninc = ( dx * 360. ) / 2.0 / np.pi / 6.37e6
    else:
        pole_lat = 90.0
        pole_lon = 0.0
        latinc = 0.0
        loninc = 0.0

    re = 6.37e6                
    rebydx = re / dx
    radperdeg = np.pi/180.0
    degperrad = 180.0/np.pi
    hemi = 1.0
    if (truelat1 < 0.0):
        hemi = -1.0

    if (map_proj == 3): # mercator
        clain = np.cos(radperdeg * truelat1)
        dlon = dx / ( re * clain)
        rsw = 0.0
        if (ref_lat != 0):
                rsw = np.log(np.tan(0.5*((ref_lat+90.0) * radperdeg))) / dlon
        deltalon = lon - ref_lon
        if (deltalon < -180.0):
                deltalon = deltalon + 360.0
        if (deltalon > 180.0):
                deltalon = deltalon - 360.0
        i = 0 + (deltalon / (dlon * degperrad))
        j = 0 + np.log(np.tan(0.5*((lat+90.0)*radperdeg)))/dlon - rsw

    elif (map_proj == 2): # polar-stereo
        reflon = stand_lon + 90.0
        scale_top = 1.0 + hemi*np.sin(truelat1*radperdeg)
        ala1 = ref_lat*radperdeg
        rsw = rebydx*np.cos(ala1)*scale_top/(1.0+hemi*sin(ala1))
        alo1 = (ref_lon - reflon)*radperdeg
        polei = 1.0 - rsw * np.cos(alo1)
        polej = 1.0 - hemi*rsw*np.sin(alo1)
        ala = lat*radperdeg
        rm = rebydx*np.cos(ala)*scale_top / (1.0+hemi*np.sin(ala))
        alo = (lon-reflon)*radperdeg
        i = polei + rm*np.cos(alo)
        j = polej + hemi*rm*np.sin(alo)

    elif (map_proj == 1): # lambert
        if (np.abs(truelat2) > 90.0):
                truelat2 = truelat1
        if (np.abs(truelat1-truelat2) > 0.1):
                cone = (np.log(np.cos(truelat1*radperdeg))-
                        np.log(np.cos(truelat2*radperdeg))) /\
                        (np.log(np.tan((90.0-np.abs(truelat1))*radperdeg*0.50))-
                         np.log(np.tan((90.0-np.abs(truelat2))*radperdeg*0.50)))
        else:
                cone = np.sin(np.abs(truelat1)*radperdeg)

        deltalon1 = ref_lon - stand_lon
        if (deltalon1 > 180.0):
                deltalon = deltalon1 - 360.0
        if (deltalon1 < -180.0):
                deltalon = deltalon1 + 360.0
        tl1r = truelat1*radperdeg
        ctl1r = np.cos(tl1r)

        rsw = rebydx*ctl1r/cone* (np.tan((90.0*hemi-ref_lat)*radperdeg/2.0) /
              np.tan((90.0*hemi-truelat1)*radperdeg/2.0))**cone

        arg = cone * (deltalon1*radperdeg)
        polei = hemi*1.0 - hemi*rsw*np.sin(arg)
        polej = hemi*1.0 + rsw*np.cos(arg)

        deltalon = lon - stand_lon
        if (deltalon > 180.0):
                deltalon = deltalon - 360.0
        if (deltalon < -180.0):
                deltalon = deltalon + 360.0

        rm = rebydx*ctl1r/cone* (np.tan((90.0*hemi-lat)*radperdeg/2.0)/
                                 np.tan((90.0*hemi-truelat1)*radperdeg/2.0))**cone
        arg = cone * (deltalon*radperdeg)
        i = polei + hemi*rm*np.sin(arg)
        j = polej - rm*np.cos(arg)

        i = hemi*i
        j = hemi*j

    #TODO: implement map projection 6
    #elif (map_proj == 6):
    else:
        print("Unsupported map projection")
        return

    i = int(i+0.5)
    j = int(j+0.5)

    return (i,j)

########################################################
