import numpy as np
from pymeteo import constants
from .. import OptionsWidget

class straight(OptionsWidget.OptionsWidget):

  def __init__(self):
    super(straight,self).__init__()
    name = 'straight (linear increase)'
    variables = [ ('z_constabv',    '6000', 'm'),
                  ('z_constblo',       '0', 'm'),
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
    zdep1 = self.getOption('z_constabv')
    zdep0 = self.getOption('z_constblo')
    umax = self.getOption('u_max')
    sf   = self.getOption('u_scaling')
    vmax = self.getOption('v_max')
    cx = self.getOption('u_adjust')
    cy = self.getOption('v_adjust')

    for k in range(len(z)):
      if (z[k] < zdep0): # constant below this height
        u[k] = 0
      elif (z[k] < zdep1): # shear
        u[k] = ((z[k]-zdep0)/(zdep1-zdep0))*umax
        u[k] = u[k]* (1 + (sf-1)*((z[k]-zdep0)/(zdep1-zdep0)))
      else: # constant section
        u[k] = umax*sf

    v[:] = vmax

    u[:] = u[:] - cx
    v[:] = v[:] - cy

    #emit
    return (z,u,v)
    
