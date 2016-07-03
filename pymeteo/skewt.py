#!/usr/bin/env python
"""
.. module:: pymeteo.skewt
   :platform: Unix, Windows
   :synopsis: Skew-T/Log-P plotting

.. moduleauthor:: Casey Webster <casey.webster@gmail.com>

This module allows plotting Skew-T/Log-P diagrams 
and hodographs from arrays of 
data with helper functions to plot directly from CM1 output files
in either Grads or HDF5 format and also from sounding data text
files in the format used for WRF and CM1 initialization.

This module also has code to produce rudimentary analysis of
the soundings data.  Currently this is limited to CAPE (surface 
and most unstable), CIN, wind shear, storm relative helicity, lifted index
and a rough estimate of storm motion based on the Bunkers (2000)
method.

.. figure:: _static/images/skewt.png
  :align: center

  Example Skew-T/Log-P with hodograph

Plotting
++++++++

High level plotting functions
-----------------------------

* :py:func:`plot_cm1h5` -- plots skewt from CM1 generated hdf5 files processed by ingest_cm1
* :py:func:`plot_cm1` -- plots skewt from CM1 generated GRaDS style output
* :py:func:`plot_wrf` -- plots skewt from WRF generated NetCDF output
* :py:func:`plot_sounding_data` -- plots skewt from CM1/WRF input sounding files
* :py:func:`plot_sounding_data_uwyo` -- plots skewt from uwyo sounding data (file or online)

These functions are called by command line scripts provided to make plotting from data files easy. 
You can invoke these through command line scripts as:

.. code-block:: bash

   # From tabular sounding data (e.g. WRF or CM1 intial sounding data)
   $ skewt tabular -f sounding.dat skewt.pdf

   # From GrADS stlye CM1 output
   $ skewt cm1 -p . -d cm1out -x 0 -y 0 skewt-cm1.pdf

   # From HDF5 CM1 output
   $ skewt cm1hdf -f model-data.h5 -x 0 -y 0 skewt.pdf

   # From WRF output 
   $ skewt wrf -f wrfout.nc --lat 30 --lon -80 -t 0 skewt.pdf

   # From University of Wyoming data file
   $ skewt uwyo -f uwyo-data.dat skewt.pdf

   # From University of Wyming website data (live)
   $ skewt uwyoweb --station 72251 skewt.pdf

* :py:func:`plot` -- generic high level plot function

.. code-block:: python

   import numpy as np
   import pymeteo.skewt as skewt

   # prepare 1D arrays height (z), pressure (p), potential temperature (th), 
   # water vapor mixing ratio (qv), winds (u and v) all of the same length.

   skewt.plot(None, z, th, p, qv, u, v, 'skewt.pdf')

Variables affecting the plot output
-----------------------------------

These variables affect the skewness of the plot and the bounds of the plot

* :py:data:`skew_angle` -- skewness
* :py:data:`Tmin`, :py:data:`Tmax` -- T dimension bounds (temperature, C)
* :py:data:`pbot`, :py:data:`ptop` -- p dimension bounds (pressure, Pa)

These variables affect the lines plotted on the skew-t

* :py:data:`isobars`
* :py:data:`isotherms`
* :py:data:`dry_adiabats`
* :py:data:`moist_adiabats`
* :py:data:`mixing_ratios`

These variables affect the plotting of the lines above

* :py:data:`lw_major`, :py:data:`lw_minor`
* :py:data:`lc_major`, :py:data:`lc_minor`

These variables affect the plotting style of skewt data

* :py:data:`linecolor_T`, :py:data:`linewidth_T`
* :py:data:`linecolor_Tve`, :py:data:`linewidth_Tve`, :py:data:`linestyle_Tve`
* :py:data:`linecolor_Td`, :py:data:`linewidth_Td`
* :py:data:`linecolor_Twb`, :py:data:`linewidth_Twb`
* :py:data:`linecolor_Parcel_T`, :py:data:`linewidth_Parcel_T`
* :py:data:`linecolor_Tvp`, :py:data:`linewidth_Tvp`, :py:data:`linestyle_Tvp`


Functions to draw isolines
--------------------------

If you are manually plotting skew-t data, these functions can be used to draw
various lines on your plot.  These are used by the high level plotting functions.

* :py:func:`draw_isotherms`
* :py:func:`draw_isobars`
* :py:func:`draw_dry_adiabat`
* :py:func:`draw_moist_adiabat`
* :py:func:`draw_water_mix_ratio`

Module reference
++++++++++++++++

"""
import os
import math
import pymeteo.thermo as met
import pymeteo.dynamics as dyn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import h5py
import pymeteo.cm1.read_grads as cm1
import pymeteo.interp
import pymeteo.constants as metconst
from netCDF4 import Dataset
import datetime as dt
import pymeteo.wrf as wrf
import datetime
import pymeteo.uwyo as uwyo

# This defines the skew-angle of the T axis
skew_angle = 37.5 
"""This defines the skewness of the T axis"""

# These define the domain plotted.  T values are at @1000 mb (e.g. unskewed)
Tmin = -40.0
"""Sets the left boundary of the plot. Temperature at 1000 mb (C)"""
Tmax = 40.0
"""Sets the right boundary of the plot. Temperature at 1000 mb (C)"""
pbot = 105000.0
"""Sets the bottom boundary of the plot. Pressure (Pa)"""
ptop = 10000.0
"""Sets the top boundary of the plot. Pressure (Pa)"""

## Values below used for plotting 
dp = 5000.0
"""The delta pressure used in calculating adiabats"""
ptickbot = 100000.0
"""Lowest elevated pressure level to be labelled in the plot (Pa)"""
pticktop = 10000.0
"""Highest elevated pressure level to be labelled in the plot (Pa)"""
tickdp = 10**4
"""The spacing between pressure labels in the plot (Pa)"""
plevs = np.arange(pbot,ptop-1,-dp)
"""List of pressure levels to do calculations on"""
fontscalefactor = 1

## Lines to plot
isotherms    = np.arange(-150,51,10)              # in degrees C
"""List of isotherms to plot on the skew-t.  In degrees C"""
isobars      = np.arange(ptickbot,ptop-1,-5000)   # in Pa
"""List of isobars to plot on the skew-t.  In Pascals"""
dry_adiabats = np.arange(-40,210,10)              # in degrees C
"""List of dry adiabats to plot on the skew-t.  In degrees C"""
moist_adiabats = np.concatenate((np.arange(-15.,10.1,5.),np.arange(12.,45.1,2.)))
"""List of moist adiabats to plot on the skew-t.  In degrees C"""
mixing_ratios = [0.2,0.4,0.8,1,2,3,4,6,8,10,14,18,24,32,40]
"""List of water vapor mixing ratio lines to plot. In g/kg"""

## Linewidths
lw_major = 0.6
"""Line width of 'major' lines.  E.g. Lines plotted at 10 C intervals or 50 mb intervals"""
lw_minor = 0.25
"""Line width of 'minor' lines.  E.g. Lines not plotted at the major intervals"""

## Linecolors
lc_major = 'grey'
"""Line color of 'major' lines.  E.g. Lines plotted at 10 C intervals or 50 mb intervals"""
lc_minor = 'lightgrey'
"""Line color of 'minor' lines.  E.g. Lines not plotted at the major intervals"""

## Skew-T line parameters
linecolor_T = 'black'
"""Line color of environmental temperature profile"""
linewidth_T = 1.5
"""Line width of environmental temperature profile"""

linecolor_Td = 'green'
"""Line color of environmental dew-point temperature profile"""
linewidth_Td = 1.5
"""Line width of environmental dew-point temperature profile"""

linecolor_Parcel_T = 'red'
"""Line color of lifted surface parcel temperature profile"""
linewidth_Parcel_T = 1.0
"""Line width of lifted surface parcel temperature profile"""

linecolor_Twb = 'blue'
"""Line color of environmental wet-bulb temperature profile"""
linewidth_Twb = 0.5
"""Line width of environmental wet-bulb temperature profile"""

linecolor_Tve = 'black'
"""Line color of environmental virtual temperature profile"""
linewidth_Tve = 0.7
"""Line width of environmental virtual temperature profile"""
linestyle_Tve = '--'
"""Line style of environmental virtual temperature profile"""

linecolor_Tvp = 'red'
"""Line color of lifted surface parcel virtual temperature profile"""
linewidth_Tvp = 0.7
"""Line width of lifted surface parcel virtual temperature profile"""
linestyle_Tvp = '--'
"""Line style of lifted surface parcel virtual temperature profile"""

#for plotted lines
pb_plot=105000
pt_plot=10000
pt_plot2=20000
dp_plot=1000
plevs_plot = np.arange(pb_plot,pt_plot-1,-dp_plot)
plevs_plot2 = np.arange(pb_plot,pt_plot2-1,-dp_plot)
plevs_std = [100000,85000,70000,50000,40000,30000,25000,20000,15000]

#TODO: enforce square w/ aspect ratio of plot
#Domain of the hodograph
umin = -22.5
umax = 27.5
vmin = -12.5
vmax = 27.5

##################################################################################
def plot_cm1h5(filename, xi, yi, output):
    """ Plots a skewt from an HDF5 file.

    :param filename: The name of the HDF5 file to open.
    :type filename: str
    :param xi: The X gridpoint of the skewt plot.
    :type xi: int
    :param yi: The Y gridpoint of the skewt plot.
    :type yi: int
    :param output: Output filename to save plot
    :type output: str

    To use this function the HDF5 file must have the following
    variables:

    - *xh* -- x scalar grid locations
    - *yh* -- y scalar grid locations
    - *z* -- z grid locations
    - *time* -- timestep of file data (scalar)
    - *th* -- potential temperature (K)
    - *thpert* -- potential temperature perturbation (K)
    - *prs* -- pressure (Pa)
    - *u* -- u wind speed (m/s)
    - *v* -- v wind speed (m/s)
    - *qv* -- water vapor mixing ratio (kg/kg)

    The names of these variables correspond to default naming by 
    CM1 using HDF5 output.

    """
    f = h5py.File(filename, 'r')
    z = f["/mesh/zh"][:]    # m 

    x = f["/mesh/xh"][xi]   # m
    y = f["/mesh/yh"][yi]   # m
    t = f["/time"][0]       # s
    th = f["/3d_s/thpert"][:,yi,xi] + f["/basestate/th0"][:] # K
    p = f["/3d_s/ppert"][:,yi,xi] + f["/basestate/pres0"][:] # Pa
    u = f["/3d_u/u"][:,yi,xi] # m/s
    v = f["/3d_v/v"][:,yi,xi] # m/s
    qv = f["/3d_s/qvpert"][:,yi,xi] + f["/basestate/qv0"][:] #kg/kg

    print(x,y,z[0],t,th[0],u[0],v[0],p[0],qv[0])
    plot_old(x,y,z,t,th,p,qv,u,v,filename, output)

##################################################################################
def plot_wrf(filename, lat, lon, time, output):
    """ Plots a skewt from an WRF NetCDF output file.

    :param filename: The name of the NetCDF file to open.
    :type filename: str
    :param lat: The latitude of the location of the skewt plot.
    :type lat: float
    :param lon: The longitude of the location of the skewt plot.
    :type lon: int
    :param output: Output filename to save plot
    :type output: str

    This function assumes the NetCDF file was produce by the WRF model.
    If another program produces the file but names variables in the same
    manner as the WRF model, this function should work for that dat as well.

    The variables using in the wrfout.nc data are:
    - Times ?
    - U(Time, bottom_top, south_north, west_east_stag)
    - V(Time, bottom_top, south_north_stag, west_east)
    - T(Time, bottom_top, south_north, west_east)   Theta perturbation
    - P(Time, bottom_top, south_north, west_east)   Pres perturbation
    - PB(Time, bottom_top, south_north, west_east)  pres base state
    - QVAPOR(Time, bottom_top, south_north, west_east)

    - PH(Time, bottom_top_stag, south_north, west_east)
    - PHB(Time, bottom_top_stag, south_north, west_east)
      Z = (((PH(k)+PH(k+1)) / 2) + ((PHB(k)+(PHB(k+1)) / 2) / 9.81

    - TH2(Time, south_north, west_east) -- 2m pot temp
    - PSFC(Time, south_north, west_east) -- surface pres
    - U10(Time, south_north, west_east) -- 10m U
    - V10(Time, south_north, west_east) -- 10m V
    - XTIME(Time) -- minutes since start of sim

    - XLAT, XLONG, XLAT_U/V, XLONG_U/V
    - Start_date

    Need time, theta, pressure
    """

    # Open NetCDF file
    f = Dataset(filename, 'r')

    wrf_lats = f.variables['XLAT'][time,:,:] # :,0
    wrf_lons = f.variables['XLONG'][time,:,:]# 0,:
    wrf_time = (f.variables['Times'][time]).tostring().decode('UTF-8')

    # refactor this into a wrf utility module

    map_proj = f.getncattr('MAP_PROJ')
    truelat1 = f.getncattr('TRUELAT1')
    truelat2 = f.getncattr('TRUELAT2')
    stand_lon= f.getncattr('STAND_LON')
    dx = f.getncattr('DX')
    dy = f.getncattr('DY')
    ref_lat = wrf_lats[0,0]
    ref_lon = wrf_lons[0,0]

    i,j = wrf.ll_to_ij(map_proj, truelat1, truelat2, stand_lon, dx, dy, ref_lat, ref_lon, lat, lon)

    ref_lat_u = f.variables['XLAT_U'][time,0,0]
    ref_lon_u = f.variables['XLONG_U'][time,0,0]
    ref_lat_v = f.variables['XLAT_V'][time,0,0]
    ref_lon_v = f.variables['XLONG_V'][time,0,0]

    i_u,j_u = wrf.ll_to_ij(map_proj, truelat1, truelat2, stand_lon, dx, dy, ref_lat_u, ref_lon_u, lat, lon)
    i_v,j_v = wrf.ll_to_ij(map_proj, truelat1, truelat2, stand_lon, dx, dy, ref_lat_v, ref_lon_v, lat, lon)

    # location
    N = 'N'
    if (lat < 0.):
        N = 'S'
    E = 'E'
    if (lon < 0.):
        E = 'W'
    x = "{0} {1}, {2} {3}".format(abs(lat), N, abs(lon), E)
    
    # pressure
    p_surface = f.variables['PSFC'][time,j,i]
    p = f.variables['P'][time,:,j,i] + f.variables['PB'][time,:,j,i]
    p = np.insert(p, 0, p_surface)

    # z heights
    ph = f.variables['PH'][time,:,j,i]
    phb = f.variables['PHB'][time,:,j,i]
    z_surface = ph[0]
    z = np.zeros(len(ph)-1, np.float32)
    for k in range(len(z)):
         z[k] = ((( ph[k]+ph[k+1] ) / 2.0) + ( phb[k] + phb[k+1]) / 2.0 ) / 9.81
    z = np.insert(z, 0, z_surface)

    # t
    t = wrf_time

    # theta
    th_surface = f.variables['TH2'][time,j,i]
    th = f.variables['T'][time,:,j,i] + 300.0
    th = np.insert(th, 0, th_surface)

    # u
    u_surface = f.variables['U10'][time,j,i]
    u = f.variables['U'][time,:,j_u,i_u]
    u = np.insert(u, 0, u_surface)

    # v
    v_surface = f.variables['V10'][time,j,i]
    v = f.variables['V'][time,:,j_v,i_v]
    v = np.insert(v, 0, v_surface)

    # qv
    qv_surface = f.variables['Q2'][time,j,i]
    qv = f.variables['QVAPOR'][time,:,j,i]
    qv = np.insert(qv, 0, qv_surface)

    title = os.path.basename(filename)
    print(x, z[0],t,th[0],u[0],v[0],p[0],qv[0])
    plot(x, z, th, p, qv, u, v, output, t, title)

    
##################################################################################
#
def plot_sounding_data(filename, output):
        """Plot SkewT from a WRF / CM1 compatible sounding data file
    
        :param filename: The name of the file to open.
        :type filename: str
        :param output: The name of the file to output plot
        :type output: str

        The datafile is the same format as used in sounding initalization
        files for WRF and CM1. 
        
        The format is:
        1 line header with surface pressure (mb), theta (K) and qv (g/kg)
        n number of lines with z (m), theta (K), qv (g/kg), u (m/s), v(m/s)

        """
        # load first line of file
        with open(filename, 'r') as f:
            surface = f.readline()

        p0, th0, qv0 = surface.split()

        # load rest of file
        _z, _th, _qv, _u, _v = np.loadtxt(filename, unpack=True, skiprows=1)

        # create arrays with one more z index for surface values
        nk = len(_z) + 1
        z = np.empty(nk, np.float32)
        th= np.empty(nk, np.float32)
        qv= np.empty(nk, np.float32)
        u = np.empty(nk, np.float32)
        v = np.empty(nk, np.float32)
        p = np.empty(nk, np.float32)

        # copy the arrays, leaving room at the surface
        z[1:nk] = _z
        th[1:nk] = _th
        qv[1:nk] = _qv / 1000.
        u[1:nk] = _u
        v[1:nk] = _v

        # assign surface values
        z[0] = 0.
        th[0] = float(th0)
        qv[0] = float(qv0) / 1000.
        u[0] = 1.75*u[1]-u[2]+0.25*u[3]
        v[0] = 1.75*v[1]-v[2]+0.25*v[3]
        p[0] = float(p0) * 100.

        # integrate pressure, assume hydrostatic
        # dp = -rho*g*dz
        for k in np.arange(1,nk):
            p[k] = p[k-1] * np.exp((metconst.g*(z[k-1]-z[k]))/(metconst.Rd*met.T((th[k]+th[k-1])/2.,p[k-1])))

        #for k in np.arange(nk):
        #    print(z[k], p[k], th[k], qv[k], u[k], v[k])
            
        plot(None, z, th, p, qv, u, v, output, title="input sounding")

def plot_sounding_data_uwyo(filename, output, stationID=0, date=None):
    """Plot SkewT from University of Wyoming sounding data

    :param filename: The name of the file containing sounding data
    :type filemane: str
    :param output: The name of the file to output plot
    :type output: str
    :param stationID: The station ID of a sounding station
    :type stationID: int
    :param date: Date and time of a sounding to request
    :type date: datetime

    If filename is not None then this function will plot data from that file.
    If filename is None and stationID is non-zero, then this will request a sounding
    datafile from http://weather.uwyo.edu for the specified date.  If date is None 
    then the requested sounding will be the most recent sounding taken at either
    12Z or 00Z. 

    """

    if (filename != None):
        title, p, z, qv, wind_dir, wind_speed, th = uwyo.fetch_from_file(filename)
    elif (stationID != 0):
        if date is None:
            date = datetime.datetime.utcnow()
        title, p, z, qv, wind_dir, wind_speed, th = uwyo.fetch_from_web(date, stationID)
    else:
        print("Neither input file or station ID was provided.  No output.")
        
    p, z, qv, u, v, th = uwyo.transform_and_check_data(p, z, qv, wind_dir, wind_speed, th)    
    plot(None, z, th, p, qv, u, v, output, title=title)
    
    
        
def plot_sounding_data_csv(filename, output):
        """Plot SkewT from a CSV sounding data file
    
        :param filename: The name of the file to open.
        :type filename: str
        :param output: The name of the file to output plot
        :type output: str

        The datafile format is CSV with the following columns: 

        - pressure (mb)
        - height (m)
        - temperature (C)
        - dew point (C)
        - wind direction (degrees)
        - wind speed (m/s)
        
        Missing values should be filled with the value -9999.00
        """
        p,z,T,Td,wdir,wspd = np.loadtxt(filename, delimiter=',',  unpack=True)
        # Pressure to Pa
        p = p * 100.

        # z to km
        #z = z / 1000.
        
        # interpolate missing wind

        nk = len(z)
        u = np.empty(nk, np.float32)
        v = np.empty(nk, np.float32)

        for k in range(nk):
           if wdir[k] == -9999. and wspd[k] == -9999.:
              u[k] = v[k] = -9999.
           else:
              u[k], v[k] = dyn.wind_deg_to_uv(wdir[k], wspd[k])

           #print('{0:5.2f} {1:5.2f} = {2:5.2f} {3:5.2f}'.format(wdir[k], wspd[k], u[k], v[k]))
           
        _z = np.empty(2,np.float32)
        _u = np.empty(2,np.float32)
        _v = np.empty(2,np.float32)
        print('INTERPOLATING')
        for k in range(nk):
           if wdir[k] == -9999. and wspd[k] == -9999.:
              kb = ke = k
              while kb >= 0 and wdir[kb] == -9999. and wspd[kb] == -9999.:
                 kb -= 1
              while ke <= nk-1 and wdir[ke] == -9999. and wspd[ke] == -9999.:
                 ke += 1

              # everything in bounds
              if kb >= 0 and ke <= nk-1:
                 _z[0] = z[kb]
                 _z[1] = z[ke]
                 _u[0] = u[kb]
                 _u[1] = u[ke]
                 _v[0] = v[kb]
                 _v[1] = v[ke]
                 
                 u[k] = pymeteo.interp.linear(_z, _u, z[k])
                 v[k] = pymeteo.interp.linear(_z, _v, z[k])

              elif kb < 0:
                 u[k] = u[ke]
                 v[k] = v[ke] 
              elif ke > nk-1:
                 u[k] = u[kb]
                 v[k] = v[kb] 

        for k in range(nk):
           # kt to m/s
           u[k] = u[k] * 0.5144444
           v[k] = v[k] * 0.5144444
        #   print('{0:5.2f} {1:5.2f} = {2:5.2f} {3:5.2f}'.format(wdir[k], wspd[k], u[k], v[k]))

        # calc theta
        th = np.empty(nk, np.float32)
        # calc qv
        qv = np.empty(nk, np.float32)
        for k in range(nk):
           th[k] = met.theta(T[k]+met.T00, p[k]) 
           w = met.es(Td[k]+met.T00) / met.es(T[k]+met.T00)
           pp = met.es(T[k]+met.T00) / p[k]
           qv[k] = 0.622 * pp * w
           #qv[k] = met.es(Td[k]+met.T00) / (met.Rv * (T[k]+met.T00))

        #print(z, th, p, qv, u, v)

        plot(None, z, th, p, qv, u, v, output, title='Sounding Data')

        
##################################################################################
# plot_cm1
#
# This plots a skewt at domain location xi,yi at time t=0 for a given CM1 dataset
# in grads format
#
def plot_cm1(path, filename, xi, yi,output):
  """Plot skewt from a Grads format CM1 output file

  :parameter _path: Path to CM1 dataset
  :type path: str
  :parameter _filename: Filename of dataset
  :type filename: str
  :parameter x1: X gridpoint to plot SkewT
  :parameter y1: Y gridpoint to plot SkewT
  :parameter output: Filename to save skewt plot

  This function plots a skewt from a CM1 output file.  
  This routine uses winds interpolated to the scalar
  gridpoints.  
  """
  f = cm1.CM1(path, filename)
  _z = f.dimZ[:] * 1000.   # km->m

  x = f.dimX[xi]
  y = f.dimY[yi]
  t = int(f.dimT[0])
        
  f.read3dMultStart(t)
  _th = f.read3dMult('th')
  _p = f.read3dMult('prs')
  _u = f.read3dMult('uinterp')
  _v = f.read3dMult('vinterp')
  _qv = f.read3dMult('qv')
  f.read3dMultStop()
  nk = len(_p[xi,yi,:])
  th = np.empty(nk+1, np.float32)
  p  = np.empty(nk+1, np.float32)
  u  = np.empty(nk+1, np.float32)
  v  = np.empty(nk+1, np.float32)
  qv = np.empty(nk+1, np.float32)
  z  = np.empty(nk+1, np.float32)
  th[1:] = _th[xi,yi,:]
  p[1:] = _p[xi,yi,:]
  u[1:] = _u[xi,yi,:]
  v[1:] = _v[xi,yi,:]
  qv[1:] = _qv[xi,yi,:]
  z[1:] = _z[:]

  # surface values
  p[0] = 100000.
  th[0]= 300.
  qv[0]= 0.014
  u[0] = u[1]
  v[0] = v[1]
  z[0] = 0.

  plot_old(x,y,z,t,th,p,qv,u,v,filename, output)

##################################################################################
# plot
#
# This is the main skewT plotting function for a single output page containing
# A skewt, hodograph and an information block (currently disabled).
#
#TODO: turn x,y into "location text" at end of arg list default None
#TODO: put title at end default None.
#TODO: pass both to plot_the_rest
def plot_old(x, y, z, time, th, p, qv, u, v, title, output):
  plot("{0} km, {1} km".format(x,y), z, th, p, qv, u, v, output, time, title) 

def plot(loc, z, th, p, qv, u, v, output, time = None, title = None):
  """Plots Skew-T/Log-P diagrapms with hodograph

  The helper functions above facilitate loading data from
  various formats and then call this function.  If you have
  data in another format or arrays of data in python already,
  then this is the function you want to use.

  :parameter loc: Location string
  :parameter z: z grid mesh (1D)
  :parameter time: Time string
  :parameter th: Potential temperature at z points
  :parameter p: Pressure at z points
  :parameter qv: Water vapor mixing ratio at z points
  :parameter u: u winds at z points
  :parameter v: v winds at z points
  :parameter title: Title for plot
  :parameter output: Filename to save plot to

  """
  fig = plt.figure(1, figsize=(10, 8), dpi=300, edgecolor='k')
  
  # sounding
  ax1 = plt.subplot(121)
  plot_sounding_axes(ax1)
  plot_sounding(ax1, z, th, p, qv, None, None)
  # hodograph
  ax2 = plt.subplot(222)
  plot_hodo_axes(ax2)
  plot_hodograph(ax2, z, u, v)
  # datablock
  ax3 = fig.add_subplot(224)
  try:
    plot_datablock(ax3, loc, z, time, th, p, qv, u, v, title)
  except:
      print("Error calcualting sounding stats, datablock omitted");
    
  # wind barbs
  ax4 = fig.add_subplot(132)
  plot_wind_axes(ax4)
  plot_wind_barbs(ax4,z,p,u,v)
  # legend
  ax5 = fig.add_subplot(4,4,15)
  plot_legend(ax5)

  
  # Adjust plot margins.
  plt.subplots_adjust(left=0.03, bottom=0.03, right=0.97, top=0.97, wspace=0.12, hspace=0.12)
  plt.savefig(output, dpi=300,bbox_inches=0)
  plt.close()

def plot_sounding_axes(axes):
  """Plots Skew-T/Log-P axes

  This will plot isotherms, isobars, dry and moist adiabats, 
  lines of constant water vapor mixing ratio, labels and 
  setup the y axes to be reversed.

  :paramter axes: The axes to draw on
  """
  draw_isotherms(axes)
  draw_isobars(axes)
  draw_dry_adiabat(axes)
  draw_moist_adiabat(axes)
  draw_water_mix_ratio(axes)
  remove_tick_labels(axes)
  axes.axis([Tmin, Tmax, pbot, ptop])
  axes.set_ylim(axes.get_ylim()[::1])

def plot_wind_axes(axes):
  # plot wind barbs
  # TODO: also do storm-relative winds
  draw_wind_line(axes)
  axes.set_axis_off()
  axes.axis([-1,1,pbot,ptop])

def plot_hodo_axes(axes):
  """Plot hodograph axes

  This will plot range arcs and labels for a hodograph plot
  """
  bounds = [-25,25,-25,25]
  axes.axis('equal')
  draw_hodograph(axes, bounds)
  remove_tick_labels(axes)
  axes.axis(bounds)


def plot_legend(axes):
  """Plot skew-t legend"""

  tT = r'Temperature'
  lT = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=linewidth_T, color=linecolor_T)

  tTd = r'Dew-point Temperature'
  lTd = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=linewidth_Td, color=linecolor_Td)

  tPT = r'Lifted Surface Parcel Temperature'
  lPT = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=linewidth_Parcel_T,
               color=linecolor_Parcel_T)

  tTwb = r'Wet-bulb Temperature'
  lTwb = Line2D(range(10), range(10), linestyle='-', marker='', linewidth=linewidth_Twb, color=linecolor_Twb)

  tTve = r'Virtual Temperature'
  lTve = Line2D(range(10), range(10), linestyle=linestyle_Tve, marker='', linewidth=linewidth_Tve,
                color=linecolor_Tve)

  tTvp = r'Lifted Surface Parcel Virtual Temperature'
  lTvp = Line2D(range(10), range(10), linestyle=linestyle_Tvp, marker='', linewidth=linewidth_Tvp,
                color=linecolor_Tvp)

  plt.legend((lT, lTve, lTd, lTwb, lPT, lTvp,),(tT, tTve, tTd, tTwb, tPT, tTvp,),
             loc=(0.125,0), fontsize=6, handlelength=10)
  # loc =, frameon=, fontsize=
  axes.set_axis_off()
  
  
def plot_wind(axes, z, p, u, v, x=0):  
  for i in np.arange(0,len(z),1):
    if (p[i] > pt_plot):
      plt.barbs(x,p[i],u[i],v[i], length=5, linewidth=.5)

  
def plot_sounding(axes, z, th, p, qv, u = None, v = None):
  """Plot sounding data

  This plots temperature, dewpoint and wind data on a Skew-T/Log-P plot.
  This will also plot derived values such as wetbulb temperature and
  label the surface based LCL, LFC and EL.

  :parameter z: height values (1D array)
  :parameter th: potential temperature at z heights (1D array)
  :parameter p: pressure at z heights (1D array)
  :parameter qv: water vapor mixing ratio at z heights (1D array)
  :parameter u: U component of wind at z heights (1D array)
  :parameter v: V component of wind at z heights (1D array)
  :paramter axes: The axes instance to draw on
  """
  # calculate Temperature and dewpoint
  T = met.T(th,p) - met.T00                          # T (C)
  Td = met.Td(p, qv) - met.T00                       # Td (C)

  # calculate wetbulb temperature
  Twb = np.empty(len(z), np.float32)                  # Twb (C)
  for zlvl in range(len(z)):
    Twb[zlvl] = met.Twb(z, p, th, qv, z[zlvl])

  # Get surface parcel CAPE and temperature / height profiles
  pcl = met.CAPE(z, p, T+met.T00, qv, 1)        # CAPE
  T_parcel = pcl['t_p'] - met.T00                      # parcel T (C)
  T_vparcel = pcl['tv_p'] - met.T00                     # parcel Tv (C)
  T_venv = met.T(pcl['thv_env'], pcl['pp']) - met.T00  # Env Tv (C)

  # plot Temperature, dewpoint, wetbulb and lifted surface parcel profiles on skew axes
  axes.semilogy(T + skew(p), p, basey=math.e, color=linecolor_T , linewidth = linewidth_T)
  axes.semilogy(Td + skew(p), p, basey=math.e, color=linecolor_Td, linewidth = linewidth_Td)
  axes.semilogy(T_parcel + skew(pcl['pp']), pcl['pp'], basey=math.e,
                color=linecolor_Parcel_T, linewidth=linewidth_Parcel_T)
  axes.semilogy(Twb + skew(p), p, basey=math.e, color=linecolor_Twb, linewidth=linewidth_Twb)

  # plot virtual temperature of environment and lifted parcel
  axes.semilogy(T_venv + skew(pcl['pp']), pcl['pp'], basey=math.e, color=linecolor_Tve,
                linewidth=linewidth_Tve, linestyle=linestyle_Tve)
  axes.semilogy(T_vparcel + skew(pcl['pp']), pcl['pp'], basey=math.e, color=linecolor_Tvp,
                linewidth=linewidth_Tvp, linestyle=linestyle_Tvp)

  # Add labels for levels based on surface parcel
  #debug print(pcl['lfcprs'], pcl['lclprs'], pcl['elprs'], pcl['ptops'])
  if (pcl['lfcprs'] > 0):
    label_m(Tmax-.5, pcl['lfcprs'], '--LFC', axes)
  if (pcl['lclprs'] > 0):
    label_m(Tmax-.5, pcl['lclprs'], '--LCL', axes)
  if (pcl['elprs'] > 0):
    label_m(Tmax-.5, pcl['elprs'], '--EL', axes)
  if (pcl['ptops'] > 0):
    label_m(Tmax-.5, pcl['ptops'], '--TOPS', axes)

  # plot labels for std heights
  for plvl in plevs_std:
    zlvl = pymeteo.interp.interp_height(z,p,plvl)
    label_m(Tmin-.5,plvl, str(int(zlvl)), axes)

  # plot wind barbs on left side of plot.  move this?  right side?
  if (u is not None and v is not None):
      #draw_wind_line(axes)
      for i in np.arange(0,len(z),2):
          if (p[i] > pt_plot):
              plt.barbs(Tmin+4,p[i],u[i],v[i], length=5, linewidth=.5)

def plot_wind_barbs(axes, z, p, u, v):
    for i in np.arange(0,len(z)):
        if (p[i] > pt_plot):
            plt.barbs(0,p[i],u[i],v[i], length=5, linewidth=.5)

              
def plot_hodograph(axes, z, u, v):
  """Plot Hodograph data

  This plots u and v winds vs height on a hodograph.

  :parameter z: height values (1D array)
  :parameter u: U component of wind at z heights (1D array)
  :parameter v: V component of wind at z heights (1D array)
  :paramter axes: The axes instance to draw on
  """
  
  # plot hodograph
  z6km = 0
  while z[z6km] <= 12000:
    z6km += 1
  axes.plot(u[0:z6km],v[0:z6km], color='black', linewidth=1.5)

  for zlvl in np.arange(0,7000,1000):
    ulvl = pymeteo.interp.linear(z,u,zlvl)
    vlvl = pymeteo.interp.linear(z,v,zlvl)
    #print('calculating winds at height {0} = ({1},{2})'.format(zlvl,ulvl,vlvl))
    label_h2(ulvl+1,vlvl-1,str(zlvl/1000), 'black', 0, axes)
    axes.plot(ulvl,vlvl, color='black', markersize=5, marker='.')

  #TODO: fix this
    try:
      ucb = dyn.storm_motion_bunkers(u,v,z)
      axes.plot(ucb[0],ucb[1],markersize=4,color='black',marker='x')
      axes.plot(ucb[2],ucb[3],markersize=4,color='black',marker='x')
    except:
        print("Error calculating sounding stats, storm motion marker not plotted");
      
def calc_sounding_stats(_z, _th, _p, _qv):
  T = met.T(_th,_p)                        # T (K)
  pcl = met.CAPE(_z, _p, T, _qv, 1)        # CAPE
  mupcl = met.CAPE(_z, _p, T, _qv, 2)      # MUCAPE
  mlpcl = met.CAPE(_z, _p, T, _qv, 3)      # MLCAPE

  return (pcl,mupcl,mlpcl)

def calc_hodograph_stats(_z, _u, _v):

  try:
      ucb = dyn.storm_motion_bunkers(_u,_v,_z)

# SRH
      srh01 = dyn.srh(_u, _v, _z, 0., 1000., ucb[0], ucb[1])
      srh03 = dyn.srh(_u, _v, _z, 0., 3000., ucb[0], ucb[1])

  except:
      print("Error calculating storm motion")
      ucb = [0.,0.,0.,0.]
      srh01 = [0.,0.]
      srh03 = [0.,0.]
      
  erh01 = dyn.srh(_u, _v, _z, 0., 1000., 0., 0.)
  erh03 = dyn.srh(_u, _v, _z, 0., 3000., 0., 0.)
 
  shear01 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 0., 1000.))
  shear03 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 0., 3000.))
  shear06 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 0., 6000.))
  shear12 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 1000., 2000.))
  shear23 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 2000., 3000.))
  shear34 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 3000., 4000.))
  shear45 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 4000., 5000.))
  shear56 = dyn.uv_to_deg(*dyn.shear(_u, _v, _z, 5000., 6000.))

  dict = { 'bunkers' : ucb,
           'srh01'   : srh01,
           'srh03'   : srh03,
           'erh01'   : erh01,
           'erh03'   : erh03,
           's01'     : shear01,
           's03'     : shear03,
           's06'     : shear06,
           's12'     : shear12,
           's23'     : shear23,
           's34'     : shear34,
           's45'     : shear45,
           's56'     : shear56
         }

  return dict


def plot_datablock(ax4, _x,_z,_t,_th,_p,_qv,_u,_v, _title):
  pcl, mupcl, mlpcl = calc_sounding_stats(_z, _th, _p, _qv)
  shear = calc_hodograph_stats(_z, _u, _v)

  brn = dyn.brn(_u, _v, _z, pcl['cape'])

  # draw datablock
  ax4.set_axis_off()
  plt.axis([-1,1,-1,1])
  plt.text(0,1,_title, verticalalignment='top', horizontalalignment='center', weight='bold', fontsize=10)
  line = ""
  if (_x != None):
    line += "Sounding at location " + str(_x) + "."
  if (_t != None):
    line += " Time = " + str(_t) + "."
  line += " " + str(len(_z)) + ' vertical levels'
  plt.text(0,.85, line, verticalalignment='center', horizontalalignment='center', fontsize=5)

  cth,cr = dyn.uv_to_deg(shear['bunkers'][0],shear['bunkers'][1])
  line = 'Estimated storm motion (supercell right mover) -> {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(cth),cr)
  plt.text(-1,.75, line, verticalalignment='center', horizontalalignment='left', fontsize=5)

  print_parcel_info('Surface Parcel', pcl, -1., .65)
  print_parcel_info('Most Unstable Parcel', mupcl, -0.3, .65)
  print_parcel_info('500 m Mixed Layer Parcel', mlpcl, 0.4, .65)

	# LCL, CCL, EL, convective temp?
	# other data?
	
  x = 0.4
  y = 0
  line = 'Hodograph'
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  x += 0.02
  y -= 0.065
  line = '0-1 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s01'][0]),shear['s01'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '0-3 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s03'][0]),shear['s03'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '0-6 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s06'][0]),shear['s06'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = 'SRH 0-1 : {0:d} m2/s2'.format(int(shear['srh01']))
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = 'SRH 0-3 : {0:d} m2/s2'.format(int(shear['srh03']))
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = 'ERH 0-1 : {0:d} m2/s2'.format(int(shear['erh01']))
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = 'ERH 0-3 : {0:d} m2/s2'.format(int(shear['erh03']))
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = 'BRN : {0:d}'.format(int(brn))
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '0-1 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s01'][0]),shear['s01'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '1-2 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s12'][0]),shear['s12'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '2-3 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s23'][0]),shear['s23'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '3-4 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s34'][0]),shear['s34'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '4-5 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s45'][0]),shear['s45'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.05
  line = '5-6 km shear {0:3d}$^\circ$ {1:3.1f} m/s'.format(int(shear['s56'][0]),shear['s56'][1])
  plt.text(x,y, line, verticalalignment='center', horizontalalignment='left', fontsize=5)

def print_3col(name, value, unit, x, y):
   plt.text(x,y, name, verticalalignment='center', horizontalalignment='left', fontsize=5)
   plt.text(x+.35,y, value, verticalalignment='center', horizontalalignment='right', fontsize=5)
   plt.text(x+.4,y, unit, verticalalignment='center', horizontalalignment='left', fontsize=5)

def print_parcel_info(title, pcl, x, y):
  plt.text(x,y, title, verticalalignment='center', horizontalalignment='left', fontsize=5)
  y -= 0.065
  x += 0.02
  print_3col('CAPE', '{0}'.format(int(pcl['cape'])), 'J kg$^{-1}$', x, y)
  y -= 0.05
  print_3col('CIN', '{0}'.format(int(pcl['cin'])), 'J kg$^{-1}$', x, y)
  y -= 0.05
  print_3col('TOPS', '{0:4.2f}'.format(float(pcl['ztops'])), 'km', x, y)
  y -= 0.05
  print_3col(r'$\theta_e$', '{0:4.1f}'.format(float(pcl['theta_e'])), 'K', x, y)
  y -= 0.05
  print_3col('LI$_{MAX}$', '{0:3.1f}'.format(float(pcl['max_li'])), 'C', x, y)
  y -= 0.05
  print_3col('LI$_{500}$', '{0:3.1f}'.format(float(pcl['li500'])), 'C', x, y)
  y -= 0.05
  print_3col('LI$_{300}$', '{0:3.1f}'.format(float(pcl['li300'])), 'C', x, y)
  y -= 0.05
  print_3col('Parcel', '{0}'.format(int(pcl['prs']/100.)), 'mb', x, y)


def remove_tick_labels(axes):
  axes.tick_params(axis='x', top='off', bottom='off', which='both')#, labelsize=0)
  axes.tick_params(axis='y', left='off', right='off', which='both')#, labelsize=0)
  for xlabel_i in axes.get_xticklabels():
    xlabel_i.set_visible(False)
    xlabel_i.set_fontsize(0.0)
  for xlabel_i in axes.get_yticklabels():
    xlabel_i.set_fontsize(0.0)
    xlabel_i.set_visible(False)

def set_fontscalefactor(newfactor):
  fontscalefactor = newfactor

def label(x, y, s, c, r, axes):
	axes.text(x,y*100,s, verticalalignment='center', horizontalalignment='center', weight='bold', fontsize=5*fontscalefactor, color=c, rotation=r)

def label_h(x, y, s, c, r, axes):
	axes.text(x,y,s, verticalalignment='center', horizontalalignment='center', weight='normal', fontsize=7*fontscalefactor, color=c, rotation=r)

def label_h2(x, y, s, c, r, axes):
	axes.text(x,y,s, verticalalignment='center', horizontalalignment='center', weight='bold', fontsize=10*fontscalefactor, color=c, rotation=r)

def label_m(x, y, s, axes):
	axes.text(x,y,s, verticalalignment='center', horizontalalignment='right', fontsize=4*fontscalefactor)

#draw hodograph u-v space
def draw_hodograph(axes, bounds):
  xmin, xmax, ymin, ymax = bounds
  gmax = max(xmax,ymax)
  gmin = min(xmin,ymin)
	# draw u-v grid
  axes.plot([gmin,gmax],[0,0], color='black', linewidth=.5)
  axes.plot([0,0],[gmin,gmax], color='black', linewidth=.5)
  for u in np.arange(xmin+1,xmax):
    if (u%5==0) and (u != 0):
      label_h(u,-1,str(u),'black',0,axes)
  for v in np.arange(ymin+1,ymax):
    if (v%5==0) and (v != 0):
      label_h(-1,v,str(v),'black',0,axes)
	# draw speed rings
  for u in np.arange(5,math.sqrt(xmax**2+ymax**2)+1,5):
    axes.plot( *hodograph_circle(u), color='grey', linewidth=.3)

# helper function for circle points
def hodograph_circle(r):
	phi = np.arange(0., 2*math.pi, 0.01)
	return r*np.cos(phi), r*np.sin(phi)

# draws a vertical line axis upon which to draw wind barbs
def draw_wind_line(axes):
	wind_line = []
	for p in plevs_plot:
		wind_line.append(0)
	axes.semilogy(wind_line, plevs_plot, color='black', linewidth=.5)
	# plot circles at certain levels?
	for p in plevs_std:
		axes.semilogy(0,p, color='black', markersize=3, marker='.')

# Puts the skew in skew-T
def skew(p):
    """Puts the skew in skew-T

    :parameter p: pressure level of the point.

    This calculates the skew of the T axis for a point to plot.  
    This assumes a logarithmic y axes and uses the variable
    :py:data:skew_angle to determine the skew.  This is the 
    magic that turns a cartesian plot into a Skew-T/Log-p plot.
    """
    return skew_angle * np.log(met.p00/p)

# Draw isotherms on skew-T / log p axes
def draw_isotherms(axes):
    """Plot isotherms on axes

    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`

    This function draws isotherms every 10 C.
    """
    for T in isotherms:
        if (T % 10 == 0):
           axes.semilogy(T + skew(plevs_plot), plevs_plot, basey=math.e, color = lc_major, linewidth= lw_major)
        else:
           axes.semilogy(T + skew(plevs_plot), plevs_plot, basey=math.e, color = lc_minor, linewidth= lw_minor)
    for T in np.arange(-40, 40, 10):
        label(T+skew(87500),875, str(T), 'red', 90.-skew_angle, axes)
    for T in np.arange(-100, -20, 10):
        label(T+skew(17500),175, str(T), 'red', 90.-skew_angle, axes)

def draw_isobars(axes):
    """Plot isobars on axes

    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`

    This function draws isobars at intervals of 2*dp.
    """
    for i in isobars:
        if (i % 5000 == 0):
            axes.plot([Tmin, Tmax], [i,i], color = lc_major, linewidth = lw_major)
        else:
            axes.plot([Tmin, Tmax], [i,i], color = lc_minor, linewidth = lw_minor)
    for i in np.arange(1000,100,-50):
        label(-10-((1000-i)*.025),i,str(i),'black',0, axes)

def draw_dry_adiabat(axes):
    """Plot dry adiabats on axes

    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`

    This function calculates dry adiabats
    and plots these lines.  Adiabats are calculated 
    every 10 K
    """
    for T in dry_adiabats:
        dry_adiabat = met.T(T+met.T00,plevs_plot) - met.T00 + skew(plevs_plot)
        if (T % 10 == 0):
            axes.semilogy(dry_adiabat, plevs_plot, basey=math.e, color = lc_major, linewidth = lw_major)
        else:
            axes.semilogy(dry_adiabat, plevs_plot, basey=math.e, color = lc_minor, linewidth = lw_minor)
            
    for T in np.arange(-20, 150, 10):
        p = (600. - 3.5*T)*100.
        x = met.T(T+met.T00,p) -met.T00 + skew(p)
        x1 = met.T(T+met.T00,p+.5*dp_plot) -met.T00 + skew(p+.5*dp_plot)
        x2 = met.T(T+met.T00,p-.5*dp_plot) -met.T00 + skew(p-.5*dp_plot)
        dx = x2-x1
        theta = math.atan2(-dp_plot,-dx) * 180/math.pi +37
        label(x,p/100,str(T),'black',theta, axes)


def draw_moist_adiabat(axes):
    """Plot moist adiabats on axes

    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`

    This function calculates moist adiabats
    and plots these lines.  Adiabats are calculated for
    values of T at 1000mb from -15 to 45 C every 5 C between
    -15 and 10 C and every 2.5 C between 12.5 and 45 C.
    """
    ps_blo = [p for p in plevs_plot if p > 100000]
    ps_blo.reverse()
    ps = [p for p in plevs_plot2 if p < 100000]
    for T in moist_adiabats:
        T_1000 = T = T + met.T00
        moist_adiabat = []
        # work backwards from 1000mb
        for p in ps_blo:
            T += met.dTdp_moist(T,p) * dp_plot
            moist_adiabat.append(T - met.T00 + skew(p))
        #reverse list order
        moist_adiabat.reverse()
        # insert 1000mb point
        T = T_1000
        moist_adiabat.append(T - met.T00)
        # work forwards from 1000mb
        for p in ps:
            T -= met.dTdp_moist(T,p) * dp_plot
            moist_adiabat.append(T - met.T00 + skew(p))
            # draw labels
            if (p == 22000):
                if (T_1000 >= met.T00 and T_1000 <= 30+met.T00):
                    label(T-met.T00+skew(p),p/100,str(int(T_1000-met.T00)),'green', 0, axes)
        if (int(T_1000 - met.T00) % 5 == 0):            
            axes.semilogy(moist_adiabat, plevs_plot2, basey=math.e, color = lc_major, linewidth = lw_major)
        else:
            axes.semilogy(moist_adiabat, plevs_plot2, basey=math.e, color = lc_minor, linewidth = lw_minor)


def draw_water_mix_ratio(axes):
    """Plot lines of constant water vapor mixing ratio on axes

    :parameter axes: The axes to draw on
    :type axes: :py:class:`matplotlib.axes`

    This function calculates isolines of constant water vapor
    mixing ratio and plots these lines.  Values of w calculated
    are given by the list variable w.
    """
    #TODO: put w and the top plevel for plotting somewhere configurable
    ps = [p for p in plevs if p>=20000 and p<=105000]
    for W in mixing_ratios:
        water_mix = []
        for p in ps:
            T = TMR(W,p/100.) 
            water_mix.append(T + skew(p))
        axes.semilogy(water_mix, ps, basey=math.e, color = 'grey', linestyle = '--', linewidth = .5)

        # Label the isoline
        T = TMR(W,1075.)
        label(T+skew(107500.), 1075, str(W), 'black', -15, axes)

def label_std_heights(axes):
   """Plot heights of standard pressure levels 

   :paramter axes: The axes to draw on
   :type axes: :py:class:`matplotlib.axes`
   """
   xpos = Tmin+1.5
   std_heights = [(1000,111),(925,2512),(850,1457),(700,3012),(500,5574),(400,7185),
                 (300,9164),(250,10363),(200,11784),(150,13608)]#,(100,16180)]
   for p,z in std_heights:
      label(xpos, p, str(z), 'black', 0, axes)

def TMR(W, p):
  # Computes temperature on mixing ratio w at pressure p.
  # TMR in C, w in g/kg dry air, p in millibars.
  # TODO: change this to something else?
  x = np.log10(W * p / (622 + W))
  TMR = 10 ** (0.0498646455 * x + 2.4082965) - 280.23475 + 38.9114 * ((10 ** (0.0915 * x) - 1.2035) ** 2)
  return TMR


