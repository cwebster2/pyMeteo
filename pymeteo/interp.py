
# interpolation function

import numpy as np
import pymeteo.constants
import bisect


def linear(dim, var, dval):
   # interpolates var to dval

   if len(dim) != len(var):
      raise Exception('Dimensions of dim and var do not match')

   n = len(dim)

   #print('Looking for {0}'.format(dval))

   # dont extrapolate (set to missing?)
   if dval < dim[0]:
      return var[0]

   if dval > dim[n-1]:
      return var[n-1]

   z0 = n-1
   while dim[z0] > dval:
      z0 = z0 - 1

   #print('lower bound found idx {0} = {1}, {2}'.format(z0, dim[z0], var[z0]))

   if dim[z0] == dval:
      #print('found it')
      return var[z0]

   z1 = 0
   while dim[z1] < dval:
      z1 = z1 + 1 

   #print('upper bound found idx {0} = {1}, {2}'.format(z1, dim[z1], var[z1]))


   zdiff = dim[z1]-dim[z0]
   zdist = (dval-dim[z0]) / zdiff

   ival = var[z0]*(1-zdist) + var[z1]*(zdist)

   #print('interpolated value = {0}'.format(ival))
   return ival 


def interp_height(z, p, plvl):
   #interpolates height to a pressure level

   nlevs = len(p)

   # check bounds
   if plvl > p[0]:
      return 0
   if plvl < p[nlevs-1]:
      return -1

   z0 = 0
   while p[z0] > plvl:
      z0 = z0 + 1 

   if p[z0] == plvl:
      return z[z0]

   z1 = nlevs-1
   while p[z1] < plvl:
      z1 = z1 - 1

   #interpolate to height.  
   # Code adapted from NSHARP 95 John Hart NSSFC KCMO. 

   zdiff = z[z1] - z[z0]
   pdiff = np.log( p[z0] / p[z1])
   pdist = np.log( p[z0] / plvl)
   height = z[z0] + (( pdist / pdiff) * zdiff)

   return height

def interp_pressure(p, z, zlvl):

   nlevs = len(p)
   #print('looking for {0}'.format(zlvl))

   z0 = nlevs-1
   while z[z0] > zlvl:
     z0 = z0 - 1

   #print('lower bound found idx {0} = {1}, {2}'.format(z0, z[z0], p[z0]))

   if z[z0] == zlvl:
      return p[z0]

   z1 = 0
   while z[z1] < zlvl:
     z1 = z1 + 1
     if (z1 >= nlevs):
        return
 
   #print('upper bound found idx {0} = {1}, {2}'.format(z1, z[z1], p[z1]))

   zdiff = z[z1] - z[z0]
   zdist = zlvl - z[z0]
   pdiff = np.log( p[z1] / p[z0] )
   pres = p[z0] * np.exp((zdist/zdiff) * pdiff)
   #print(pres)
   return pres

