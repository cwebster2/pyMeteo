import numpy as np
from pymeteo import constants
from .. import OptionsWidget

class straightllws(OptionsWidget.OptionsWidget):

  def __init__(self):
    super(straightllws,self).__init__()
    name = 'straight (linear increase)'
    variables = [ ('z_constabv',    '6000', 'm'),
                  ('z_llws,top',     '1000', 'm'),
                  ('u_llws',        '12.0', 'm/s'),
                  ('u_max',         '30.0', 'm/s'),
                  ('u_scaling',      '1.0', ''),
                  ('v_max',          '7.0', 'm/s'),
                  ('u_adjust'     ,   '0.0', 'm/s'),
                  ('v_adjust'      ,  '0.0', 'm/s') ]
    self.initUI(name, variables)
 
  def plot(self):
    # nned z, t, th, p and qv
    z = np.arange(0., 22000., 50.)
    u = np.zeros(len(z))
    v = np.zeros(len(z))

    # parameters
    zdep2 = self.getOption('z_constabv')
    zdep1 = self.getOption('z_llws,top')
    ullws = self.getOption('u_llws')
    umax = self.getOption('u_max')
    sf   = self.getOption('u_scaling')
    vmax = self.getOption('v_max')
    cx = self.getOption('u_adjust')
    cy = self.getOption('v_adjust')

    for k in range(len(z)):
      if (z[k] < zdep1): # llws zone
        u[k] = (z[k]/zdep1)*ullws
      elif (z[k] < zdep2): # straight shear
        u[k] = (ullws+(z[k]-zdep1)*(umax-ullws)/(zdep2-zdep1))
        #u[k] = u[k]* (1 + (sf-1)*((z[k]-zdep0)/(zdep1-zdep0)))
      else: # constant section
        u[k] = umax

    v[:] = vmax

    u[:] = u[:] - cx
    v[:] = v[:] - cy

    #emit
    return (z,u,v)
    
