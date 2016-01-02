#!/usr/bin/env python
"""
.. module:: pymeteo.thermo
   :platform: Unix, Windows
   :synopsis: Thermodynamic routines

.. moduleauthor:: Casey Webster <casey.webster@gmail.com>

"""
#todo
# complete constants and function library

#constants

import math
import numpy as np
from pymeteo.constants import *
import pymeteo.interp

def T(theta,p):
    """Convert Potential Temperature :math:`\\theta` to Temperature
    
    :parameter theta: Potential temperature (K)
    :parameter p: Pressure (Pa)
    :returns: Temperature (K) 
    """
    return theta * (p00/p)**-kappa_d

def theta(T,p):
	return T * (p00/p)**kappa_d

def es(T):
	#return 611.2 * np.exp((Lv(T)/Rv)*((1./T00)-(1./T)))
        return 611.2 * np.exp(17.67*(T-T00)/(T-29.65))

def esi(T):
	#return 611.2 * np.exp((Lv(T)/Rv)*((1./T00)-(1./T)))
        return 611.2 * np.exp(21.8745584*(T-T00)/(T-7.66))

def w_vs(T,pd):
	return epsilon * (es(T)/pd)

def theta_v(th, qv):
   return th * (1. + 0.61*qv)

def Td(p, qv):
    old_settings = np.seterr(all='ignore')
    el = np.log((qv/epsilon)*p/100./(1. + (qv/epsilon)))
    Td = T00 + (243.5*el-440.8)/(19.48-el)
    np.seterr(**old_settings)
    return Td

def Lv(T):
#TODO: Temp dependance
	return L

def th_e(p, t, td, qv):
   if (td-t) >= -0.1:
      tlcl = t
   else:
      tlcl = 56.0 + ((td-56.0)**(-1) + 0.00125*math.log(t/td) )**(-1)

   th_e = t * ((100000./p) ** (0.2854*(1.0-0.28*qv)))*math.exp(((3376./tlcl)-2.54)*qv*(1.0+0.81*qv))

   return th_e

def q_vl(p, t):
   _es = es(t)
   q_vl = epsilon*_es/(p-_es)
   return q_vl

def q_vi(p, t):
   _es = esi(t)
   q_vi = epsilon*_es/(p-_es)
   return q_vi

# 16 Nov 2011 class notes -- PSU Meteo 531 Jerry Harrington
def dTdz_moist(T,p):
  pd = p - es(T)
  num = 1. + ((Lv(T) * w_vs(T,pd))/(Rd*T))
  den = 1. + ((Lv(T)**2 * w_vs(T,pd))/(cpd * Rv * T**2))
  return (-g/cpd)*(num/den)

def dTdp_moist(T,p):
	return dTdz_moist(T,p) * -((Rd*T)/(p*g))

def Twb(z,p,th,qv,z0):
   # returns wet bulb (deg C) for parcel from height z0

   if (len(z) != len(p) != len(_th) != len(qv)):
      raise Exception('Bounds of z, T, Td do not match')

   th_0 = pymeteo.interp.linear(z,th,z0)
   qv0 = pymeteo.interp.linear(z,qv,z0)
   p0 = pymeteo.interp.interp_pressure(p,z,z0)
   t0 = T(th_0, p0)
   td0 = Td(p0, qv0) 
   es0 = es(td0)

   #print(qv0, p0, t0, td0)
   if np.isnan(td0):
       return np.nan
   
   residual = 0.0005
   not_converged = 1
   i = 0
   Twb = td0-T00
   sign = 1
   prevsign = 1
   incr = 10
   while not_converged:
      Ewguess = es(Twb+T00)
      Eguess = Ewguess/100. - p0/100. * (t0-T00-Twb) * 0.00066 * (1. + (0.00115 * (Twb)))
      delta_e = es0/100. - Eguess

      if abs(delta_e) <= residual:
         not_converged = 0
         break

      if delta_e < 0.:
         sign = -1
      else:
         sign = 1
      
      if sign != prevsign:
         prevsign = sign
         incr /= 10.
         
#      print Twb, delta_e, incr, sign
      Twb += incr * float(sign)
      i+=1

      if i > 100:
         print('Failed to converge')
         Twb = 0.
         break

   return Twb 

def CAPE(z, p, t, q, parcel):

   #TODO: TOPS, LImax|500|300

   # T in K
   # P in Pa
   # q in kg/kg

   if (len(z) != len(p) != len(t) != len(q)):
      raise Exception('Bounds of z, T, Td do not match')

   ml_depth = .500  # for option of mixed layer parcel.
   pinc = 100. # Pa

   adiabat = 1

   debuglevel = 0

   nk = len(z)

   # goal: pi, td, th, thv 
   # have: z, p, t, q
   pi  = np.empty(nk, np.float32)
   td  = np.empty(nk, np.float32)
   th  = np.empty(nk, np.float32)
   thv = np.empty(nk, np.float32)
  
   pi = (p*rp00)**rddcp 
   td = Td(p,q)
   th = t/pi
   thv = th * (1. + reps * q)/(1. + q)

   #TODO: calc z from hydrostatic instead of using provided levels?

   # source parcel

   if (parcel == 1):
      kmax = 0

   elif (parcel == 2):
      # use most unstable parcel 

      if (p[0] < 50000.):
         kmax = 0
         maxthe = th_e(p[0], t[0], td[0], q[0])
      else:
         maxthe = 0.
         for k in range(nk):
            if p[k] < 50000.:
               break
            the = th_e(p[k], t[k], td[k], q[k])
            if (the > maxthe):
               maxthe = the
               kmax = k
      if (debuglevel >= 100):
         print('  kmax,maxthe = {0}, {1}'.format(kmax,maxthe))

   elif (parcel == 3):
      # mixed layer

      if ((z[1]-z[0]) > ml_depth):
         # second level is above mixed layer
         avgth = th[0]
         avgqv = q[0]
         kmax = 0

      elif (z[nk-1] < ml_depth):
         # all levels are in the mixed layer
         #print nk-1, z[nk-1], ml_depth
         avgth = th[nk-1]
         avgqv = q[nk-1]
         kmax = nk-1

      else:
         # calculate mixed layer properties

         avgth = 0.
         avgqv = 0.
         k = 1
         if (debuglevel >= 100):
            print('  ml_depth = {0}'.format(ml_depth))
            print('  k,z,th,q = {0}, {1}, {2}, {3}'.format(0, z[0], th[0], q[0]))

         while ((z[k] <= ml_depth) and (k < nk)):
            if (debuglevel >= 100):
               print('  k,z,th,q = {0}, {1}, {2}, {3}'.format(k, z[k], th[k], q[k]))

            avgth += 0.5 * (z[k] - z[k-1])*(th[k] + th[k-1])
            avgqv += 0.5 * (z[k] - z[k-1])*(q[k] + q[k-1])
            
            k += 1

         th2 = th[k-1]+(th[k]-th[k-1])*(ml_depth-z[k-1])/(z[k]-z[k-1])
         qv2 = q[k-1]+(q[k]-q[k-1])*(ml_depth-z[k-1])/(z[k]-z[k-1])

         if (debuglevel >= 100):
               print('  k,z,th,q = {0}, {1}, {2}, {3}'.format(999, ml_depth, th2, qv2))
            
         avgth += 0.5*(ml_depth-z[k-1])*(th2+th[k-1])
         avgqv += 0.5*(ml_depth-z[k-1])*(qv2+q[k-1])
            

         if (debuglevel >= 100):
               print('  k,z,th,q = {0}, {1}, {2}, {3}'.format(k, z[k], th[k], q[k]))
            
         avgth /= ml_depth
         avgqv /= ml_depth

         kmax = 0
 
      if (debuglevel >= 100):
         print('  avgth, avgqv = {0}, {1}'.format(avgth, avgqv))

   # define parcel properties
   narea = 0.

   if (parcel == 1) or (parcel ==2) :
      k = kmax
      th2 = th[kmax]
      pi2 = pi[kmax]
      p2  = p[kmax]
      t2  = t[kmax]
      thv2= thv[kmax]
      qv2 = q[kmax]
      b2 = 0.
   elif (parcel == 3):
      k = kmax
      th2 = avgth
      qv2 = avgqv
      thv2= th2*(1.+reps*qv2)/(1.+qv2)
      pi2 = pi[kmax]
      p2  = p[kmax]
      t2  = th2*pi2
      b2  = g*(thv2-thv[kmax])/thv[kmax]

   p0 = p2
   t0 = t2
   qv0 = qv2

   ql2 = 0.
   qi2 = 0.
   qt = qv2

   cape = 0.
   cin = 0.
   lfc = 0.
   
   doit = True
   cloud = False
   if (adiabat == 1) or (adiabat == 2):
      ice = False
   else:
      ice = True

   the = th_e(p2, t2, t2, qv2)
   th_p_e = the
   if (debuglevel >= 100):
      print('th_e = {0}'.format(the))

   pt = np.empty(nk, np.float32)
   pb = np.empty(nk, np.float32)
   pc = np.empty(nk, np.float32)
   pn = np.empty(nk, np.float32)
   ptv = np.empty(nk, np.float32)
   ptd = np.empty(nk, np.float32)
   pqv = np.empty(nk, np.float32)
   pql = np.empty(nk, np.float32)
   
   pt[k] = t2
   if cloud:
      ptd[k] = t2
   else:
      ptd[k] = Td(p2, qv2)
   ptv[k] = t2*(1.+reps*qv2)/(1.+qv2)
   pb[k] = 0.
   pqv[k] = qv2
   pql[k] = 0.

   zlcl = -1.
   zlfc = -1.
   zel = -1.
   ztops = -1.
   ptops = 0.
   plcl = 0.
   plfc = 0.
   pel = 0.
   max_li = 40.

   # Parcel Ascent starts here!
   
   if (debuglevel >= 100):
      print('  Start loop:')
      print('  p2,th2,qv2 = {0}, {1}, {2}'.format(p2, th2, qv2))

   while ((doit) and (k < (nk-1))):
      k = k+1
      b1 = b2

      dp = p[k-1]-p[k]
      nloop = 1
      if (dp >= pinc):
         nloop = 1 + int(dp/pinc)
         dp = dp/float(nloop)

      for n in range(nloop):

         p1 = p2
         t1 = t2
         pi1 = pi2
         th1 = th2
         qv1 = qv2
         ql1 = ql2
         qi1 = qi2
         thv1 = thv2

         p2 = p2 - dp
         pi2 = (p2*rp00)**rddcp

         thlast = th1
         i = 0
         not_converged = True

         while (not_converged):

            i += 1
            t2 = thlast*pi2
            if ice:
               fliq = max(min((t2-233.15)/(T00-233.15),1.),0.)
               fice = 1.-fliq
            else:
               fliq = 1.
               fice = 0.
            
            qv2 = min( qt, fliq*q_vl(p2,t2) + fice*q_vi(p2,t2))
            qi2 = max( fice*(qt-qv2), 0.)
            ql2 = max( qt-qv2-qi2, 0.)

            tbar = 0.5*(t1+t2)
            qvbar = 0.5*(qv1+qv2)
            qlbar = 0.5*(ql1+ql2)
            qibar = 0.5*(qi1+qi2)

            lhv = lv1-lv2*tbar
            lhs = ls1-ls2*tbar
            lhf = lhs - lhv

            rm = Rd + Rv*qvbar
            cpm = cp + cpv*qvbar + cpl*qlbar + cpi*qibar
            th2 = th1*math.exp( lhv*(ql2-ql1)/(cpm*tbar) + lhs*(qi2-qi1)/(cpm*tbar) + (rm/cpm - Rd/cp)*math.log(p2/p1))

            if i > 90:
               print(i,th2,thlast,th2-thlast)

            if i > 100:
               print('Error: lack of convergence, stopping iteration')
               quit()

            if abs(th2-thlast) > converge:
               thlast += 0.3*(th2-thlast)
            else:
               not_converged = False

         # end not_converged
 
         # pressure increment complete
         if (ql2 >= 1.0e-10):
            cloud = True

         if cloud and (zlcl < 0.):
            zlcl = z[k-1] + (z[k]-z[k-1])*float(n)/float(nloop)
            plcl = p[k-1] + (p[k]-p[k-1])*float(n)/float(nloop)

         if (adiabat == 1) or (adiabat == 3):
            #pseudoadiabat
            qt = qv2
            ql2 = 0.
            qi2 = 0.
         elif (adiabat <= 0) or (adiabat >= 5):
            print('Unknown adiabat')
            quit()

      # end nloop

      thv2 = th2*(1.+reps*qv2)/(1.+qv2+ql2+qi2)
      b2 = g * (thv2 - thv[k])/thv[k]
      dz = -cpdg*0.5*(thv[k]+thv[k-1])*(pi[k]-pi[k-1])

      if (thv[k] - thv2) < max_li:
         max_li = thv[k] - thv2

      if ( zlcl > 0.) and (zlfc < 0.) and (b2 > 0.):
         if (b1 > 0.):
           zlfc = zlcl
           plfc = plcl
         else:
           zlfc = z[k-1]+(z[k]-z[k-1])*(0.-b1)/(b2-b1)
           plfc = p[k-1]+(p[k]-p[k-1])*(0.-b1)/(b2-b1)

      if ( zlfc > 0.) and (zel < 0.) and (b2 < 0.):
        zel = z[k-1]+(z[k]-z[k-1])*(0.-b1)/(b2-b1)
        pel = p[k-1]+(p[k]-p[k-1])*(0.-b1)/(b2-b1)
         

      the = th_e(p2, t2, t2, qv2)

      pt[k] = t2
      if (cloud):
         ptd[k] = t2
      else:
         ptd[k] = Td(p2,qv2)
      ptv[k] = t2*(1.+reps*qv2)/(1.+qv2)
      pb[k] = b2
      pqv[k] = qv2
      pql[k] = ql2

      # Cape and cin contributions
      if (b2 >= 0.) and (b1 < 0.):
         # First positive area
         ps = p[k-1]+(p[k]-p[k-1])*(0.-b1)/(b2-b1)
         frac = b2/(b2-b1)
         parea = 0.5*b2*dz*frac
         narea -= 0.5*b1*dz*(1.-frac)
         if (debuglevel >= 200):
            print('        b1,b2 = {0},{1}'.format(b1,b2))
            print('        p1,ps,p2 = {0},{1},{2}'.format(p[k-1],ps,p[k]))
            print('        frac = {0}'.format(frac))
            print('        parea = {0}'.format(parea))
            print('        narea = {0}'.format(narea))

         cin += narea
         narea = 0.0
      elif (b2 < 0.) and (b1 > 0.):
         # first negative area
         ps = p[k-1]+(p[k]-p[k-1])*(0.-b1)/(b2-b1)
         frac = b1/(b1-b2)
         parea = 0.5*b1*dz*frac
         narea = -0.5*b2*dz*(1.-frac)
         if (debuglevel >= 200):
            print('        b1,b2 = {0},{1}'.format(b1,b2))
            print('        p1,ps,p2 = {0},{1},{2}'.format(p[k-1],ps,p[k]))
            print('        frac = {0}'.format(frac))
            print('        parea = {0}'.format(parea))
            print('        narea = {0}'.format(narea))

      elif (b2 < 0.):
         #still neg
         parea = 0.
         narea -= 0.5*dz*(b1+b2)
      else:
         parea = 0.5*dz*(b1+b2)
         narea = 0.

      cape += max(0., parea)
      pc[k] = cape
      #print parea,narea,cape,cin,zel,ztops,z[k],(b1-b2),(narea-cape)
      if (debuglevel >= 200):
         print(p2, b1, b2, cape, cin, cloud)

      if (narea >= cape) and (zel > 0.) and (ztops < 0.):
         #print('Tops')
         ztops = z[k-1]+(z[k]-z[k-1])*(b1-b2)/(narea-cape)
         ptops = p[k-1]+(p[k]-p[k-1])*(b1-b2)/(narea-cape)

      if (p[k] <= 10000.) and (b2 < 0.):
         doit = False
   
   # end parcel ascent loop
   
   li500 = 0
   li300 = 0

   lev = 50000.
   k = pymeteo.interp.interp_height(z, p, lev)
   t_p_lev = pymeteo.interp.linear(z, ptv, k)
   t_e_lev = T(pymeteo.interp.linear(z, thv, k),lev)
   li500 = t_e_lev - t_p_lev
   
   lev = 30000
   k = pymeteo.interp.interp_height(z, p, lev)
   t_p_lev = pymeteo.interp.linear(z, ptv, k)
   t_e_lev = T(pymeteo.interp.linear(z, thv, k),lev)
   li300 = t_e_lev - t_p_lev

   #print('CAPE = {0}'.format(cape))
   dict = { 'lfc' : zlfc,
            'lcl' : zlcl,
            'el'  : zel,
            'lfcprs' : plfc,
            'lclprs' : plcl,
            'elprs'  : pel,
            'ztops'  : ztops,
            'ptops'  : ptops,
            'lfcth' : 0,
            'lclth' : 0,
            'elth'  : 0,
            'zlevs' : z,
            't_p'  : pt, 
            'tv_p'  : ptv, 
            'thv_env'  : thv, 
            'pp'    : p,
            'theta_e': th_p_e, 
            'cape'  : cape,
            'cin'   : cin,
            'max_li': max_li,
            'li500' : li500,
            'li300' : li300,
            'prs'   : p[kmax]
          }

   return dict      

