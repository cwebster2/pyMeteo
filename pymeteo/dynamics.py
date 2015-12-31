#!/usr/bin/env python
"""This module provides routine thermodynamic functions

"""
#todo
# complete constants and function library

#constants

import math
import numpy as np
from pymeteo.constants import *
import pymeteo.interp

# Mean wind in the layer [zmin,zmax]
def avg_wind(_u, _v, _z, zmin, zmax):
#TODO: np.average !
#REPLACED BY mean_wind below
        u = 0
        v = 0
        n = 0
        for i in np.arange(0,len(_z),1):
                if (_z[i] > zmin and _z[i] < zmax):
                        u += _u[i]
                        v += _v[i]
                        n += 1
        if (n == 0):
            return 0
        return (u/n, v/n)

def storm_motion_rasmussen(_u,_v,_z):
        #rasmussen (1984)
        # calculate storm motion
        u_0_500 = avg_wind(_u,_v,_z, 0., 500.)
        u_4km = avg_wind(_u,_v,_z, 3500.,4500.)

        dist_60pct = math.sqrt((u_4km[0]-u_0_500[0])**2 + (u_4km[1]-u_0_500[1])**2) * 0.6
        du = u_4km[0]-u_0_500[0]
        dv = u_4km[1]-u_0_500[1]
        theta = math.atan(dv/du)
        u60 = dist_60pct * math.cos(theta) + u_0_500[0]
        v60 = dist_60pct * math.sin(theta) + u_0_500[1]
        theta -= math.pi/2.
        dist = 8.7
        u_cr = dist * math.cos(theta) + u60
        v_cr = dist * math.sin(theta) + v60
        theta += math.pi
        u_cl = dist * math.cos(theta) + u60
        v_cl = dist * math.sin(theta) + v60

        return u_cr,v_cr,u_cl,v_cl

def storm_motion_bunkers(_u,_v,_z):
        #bunkers (2000)
        # calculate storm motion
        u_0_500 = avg_wind(_u,_v,_z, 0., 500.)
        u_0_6km = avg_wind(_u,_v,_z, 0., 6000.)
        du = u_0_6km[0] - u_0_500[0]
        dv = u_0_6km[1] - u_0_500[1]
        theta = math.atan(dv/du) - math.pi/2.
        dist = 7.5
        u_cr = dist * math.cos(theta) + u_0_6km[0]
        v_cr = dist * math.sin(theta) + u_0_6km[1]
        theta += math.pi
        u_cl = dist * math.cos(theta) + u_0_6km[0]
        v_cl = dist * math.sin(theta) + u_0_6km[1]

        return u_cr,v_cr,u_cl,v_cl

def circulation(u, v, w, x, y, z):

  n = len(x)

	# local variables
  uavg = np.zeros((n), np.float32)
  vavg = np.zeros((n), np.float32)
  wavg = np.zeros((n), np.float32)

  # initialize returned variables
  vdotdl1 = np.empty((n), np.float32)
  vdotdl2 = np.empty((n), np.float32)
  C = 0.0

  if (u == missingval).any() or \
     (v == missingval).any() or \
     (w == missingval).any() or \
     (x == missingval).any() or \
     (y == missingval).any() or \
     (z == missingval).any():

     del uavg, vavg, wavg, dx, dy, dz, vdotdl1, vdotdl2
     return (missingval, 0, 0)

        # calculate wind and dl around the circuit
  uavg[0:n-1] = 0.5 * (u[0:n-1] + u[1:n])
  vavg[0:n-1] = 0.5 * (v[0:n-1] + v[1:n])
  wavg[0:n-1] = 0.5 * (w[0:n-1] + w[1:n])

  uavg[n-1] = 0.5*(u[0] + u[n-1])
  vavg[n-1] = 0.5*(v[0] + v[n-1])
  wavg[n-1] = 0.5*(w[0] + w[n-1])

        # ediff1d returns an array of the difference between
	# each element of the passed array, which is to say,
	# it returns the delta of each array element.  perfect.
  dx = np.ediff1d(x, to_end=x[0]-x[n-1])
  dy = np.ediff1d(y, to_end=y[0]-y[n-1])
  dz = np.ediff1d(z, to_end=z[0]-z[n-1])

        # assumes clockwise parcels
  C = np.sum(-uavg*dx - vavg*dy - wavg*dz) 
  vdotdl1 = - uavg*dx - vavg*dy - wavg*dz
  vdotdl2 = vdotdl1 / (dx**2 + dy**2 + dz**2)**0.5
  vdotdl1 = vdotdl1 * km2m

        
  C = C * km2m  # C now in units m2/s
        
  del uavg, vavg, wavg, dx, dy, dz
  return (C, vdotdl1, vdotdl2)


def integral_Bdz(th, thp, z):

  n = len(z)

  th_avg = np.empty((n), np.float32)
  thp_avg = np.empty((n), np.float32)
  dz = np.empty((n), np.float32)

  intBdz = 0.0

  if (th == missingval).any() or \
     (thp == missingval).any() or \
     (z == missingval).any():

     intBdz = missingval
     del th_avg, thp_avg, dz
     return intBdz

  th_avg[0:n-1] = 0.5 * (th[0:n-1] + th[1:n])
  thp_avg[0:n-1] = 0.5 * (thp[0:n-1] + thp[1:n])

  th_avg[n-1] = 0.5*(th[0] + th[n-1])
  thp_avg[n-1] = 0.5*(thp[0] + thp[n-1])

  dz = np.ediff1d(z, to_end=z[0]-z[n-1])

  intBdz = np.sum(-dz * gravity * thp_avg/(th_avg - thp_avg))

  del th_avg, thp_avg, dz
  intBdz = intBdz * km2m  # units now m2/s2

  return intBdz

def integral_dt(i, t):
   n = len(t)

   iavg = np.empty((n-1), np.float32)
   dt = np.empty((n-1), np.float32)

   iavg[0:n-1] = 0.5 * (i[0:n-1] + i[1:n])

   dt = np.ediff1d(t)

   integral = np.sum(iavg * dt)

   return integral


# helper functions

def uv_to_deg(u,v):
   """transforms u, v, to direction, maginutide

   :param u: u wind component
   :param v: v wind component
   :returns: wind direction and magnitude
   """
   direction = np.arctan2(u,v)*(180./np.pi)
   speed = (u**2 + v**2)**(0.5)

   return direction,speed

def wind_deg_to_uv(direction, speed):
   """Converts direction and speed into u,v wind

   :param direction: wind direction (mathmatical angle)
   :param speed: wind magnitude
   :returns: u and v wind components
   """
   u = speed * np.sin(np.pi * (direction+180.) / 180.)
   v = speed * np.cos(np.pi * (direction+180.) / 180.)

   return u,v

def shear(_u, _v, _z, zbot, ztop):
    """Calculates the shear in the layer between zbot and ztop

    :param _u: U winds (1D vector in z)
    :param _u: V winds (1D vector in z)
    :param _z: z heights (1D vector in z)
    :param zbot: Bottom of the layer
    :param ztop: Top of the layer

    """

    if zbot < _z[0]:
        zbot = _z[0]

    ubot = pymeteo.interp.linear(_z,_u, zbot)
    vbot = pymeteo.interp.linear(_z,_v, zbot)
    utop = pymeteo.interp.linear(_z,_u, ztop)
    vtop = pymeteo.interp.linear(_z,_v, ztop)

    u = utop-ubot
    v = vtop-vbot
    return u,v

def srh(_u,_v,_z, zbot, ztop, cx, cy):
    """Calculates the storm relative helicity in the layer between zbot and ztop

    :param _u: U winds (1D vector in z)
    :param _u: V winds (1D vector in z)
    :param _z: z heights (1D vector in z)
    :param zbot: Bottom of the layer
    :param ztop: Top of the layer
    :param cx: u component of storm motion
    :param cy: v component of storm motion

    """

    if zbot < _z[0]:
        zbot = _z[0]

    dz = 10.
    z = np.arange(zbot, ztop+dz, dz)
    nk = len(z)
    u = np.empty(nk, np.float32)
    v = np.empty(nk, np.float32)
    du = np.empty(nk-1, np.float32)
    dv = np.empty(nk-1, np.float32)
    uavg = np.empty(nk-1, np.float32)
    vavg = np.empty(nk-1, np.float32)

    for k in range(nk):
        u[k] = pymeteo.interp.linear(_z,_u,z[k]) 
        v[k] = pymeteo.interp.linear(_z,_v,z[k]) 

    du = np.ediff1d(u)
    dv = np.ediff1d(v)

    uavg[0:nk-1] = 0.5 * (u[0:nk-1] + u[1:nk])
    vavg[0:nk-1] = 0.5 * (v[0:nk-1] + v[1:nk])

    srh = np.sum(-(uavg-cx)*dv + (vavg-cy)*du)
    return srh

def mean_wind(_u,_v,_z, zbot, ztop):
    """Calculates the mean wind in the layer between zbot and ztop

    :param _u: U winds (1D vector in z)
    :param _u: V winds (1D vector in z)
    :param _z: z heights (1D vector in z)
    :param zbot: Bottom of the layer
    :param ztop: Top of the layer

    """

    if zbot < _z[0]:
        zbot = _z[0]

    dz = 10.
    z = np.arange(zbot, ztop+dz, dz)
    nk = len(z)
    u = np.empty(nk, np.float32)
    v = np.empty(nk, np.float32)

    for k in range(nk):
        u[k] = pymeteo.interp.linear(_z,_u,z[k])
        v[k] = pymeteo.interp.linear(_z,_v,z[k])

    uavg = np.mean(u, dtype=np.float64)
    vavg = np.mean(v, dtype=np.float64)
    return uavg, vavg

def brn(_u,_v,_z,cape):

   u06avg = mean_wind(_u,_v,_z,0.,6000.)
   u0500avg = mean_wind(_u,_v,_z,0.,500.)
   u = u06avg[0] - u0500avg[0]
   v = u06avg[1] - u0500avg[1]

   brn = cape / (0.5 * (u**2 + v**2))
   return brn
