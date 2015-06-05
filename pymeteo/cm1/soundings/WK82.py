import numpy as np
from pymeteo import constants
from pymeteo import thermo
from .. import OptionsWidget

class WK82(OptionsWidget.OptionsWidget):

  def __init__(self):
    super(WK82,self).__init__()
    name = 'Weisman-Klemp 82 Analytic Sounding'
    variables = [ ('z_tr',   '12000', 'm'),
                  ('theta_tr', '343', 'K'),
                  ('T_tr',     '213', 'K'),
                  ('th_0',     '300', 'K'),
                  ('qv_0',   '0.014', 'kg/kg'),
                  ('p_0',   '100000', 'Pa') ]
    self.initUI(name, variables)
 
  def plot(self):
    # nned z, t, th, p and qv
    z = np.arange(0., 22000., 50.)
    th = np.zeros(len(z))
    thv = np.zeros(len(z))
    p = np.zeros(len(z))
    qv = np.zeros(len(z))
    rh = np.zeros(len(z))
    pi = np.zeros(len(z))

    # parameters
    th_sfc = self.getOption('th_0')
    qv_pbl = self.getOption('qv_0')
    z_tr = self.getOption('z_tr')
    t_tr = self.getOption('T_tr')
    th_tr = self.getOption('theta_tr')
    p_sfc = self.getOption('p_0')

    pi_sfc = (p_sfc/constants.p00)**(constants.Rd/constants.cp)
    qv_sfc = thermo.q_vl(p_sfc,th_sfc*pi_sfc)
    thv_sfc = th_sfc*(1.+qv_sfc*constants.reps)/(1.+qv_sfc)

    rh[:] = 0.
    #calculate 
    # WK82 pp 506
    for k in range(len(z)):
      if (z[k] <= z_tr):
        th[k] = th_sfc + (th_tr - th_sfc)*((z[k]/z_tr)**(1.25))
        rh[k] = 1.0-0.75*((z[k]/z_tr)**1.25)
      else:
        th[k] = th_tr * np.exp((constants.gravity/(constants.cp*t_tr))*(z[k]-z_tr))
        rh[k] = 0.25

    qv[:] = 0.

    for n in range(20):
      for k in range(len(z)):
        thv[k] = th[k]*(1.+constants.reps*qv[k])/(1.+qv[k])
      pi[0] = pi_sfc-constants.gravity*z[0]/(constants.cp*0.5*(thv_sfc+thv[0]))
      for k in range(len(z)):
        if k == 0:
          continue
        pi[k] = pi[k-1]-constants.gravity*(z[k]-z[k-1])/(constants.cp*0.5*(thv[k]+thv[k-1]))
      for k in range(len(z)):
        p[k] = constants.p00*(pi[k])**(constants.cp/constants.Rd)
      #p[:] = constants.p00*(pi[:])**(constants.cp/constants.Rd)
      for k in range(len(z)):
        qv[k] = rh[k]*thermo.q_vl(p[k],th[k]*pi[k])
        if (qv[k] > qv_pbl):
          qv[k] = qv_pbl

    #emit
    return (z,th,p,qv)
