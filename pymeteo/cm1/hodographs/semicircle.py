import numpy as np
from pymeteo import constants
from .. import OptionsWidget

class semicircle(OptionsWidget.OptionsWidget):

  def __init__(self):
    super(semicircle,self).__init__()
    name = 'semicircle'
    variables = [ ('z_curve,top',    '6000', 'm'),
                  ('u_max',         '15.0', 'm/s'),
                  ('v_max',          '15.0', 'm/s'),
                  ('u_adjust'     ,   '0.0', 'm/s'),
                  ('v_adjust'      ,  '0.0', 'm/s') ]
    self.initUI(name, variables)
 
  def plot(self):
    # nned z, t, th, p and qv
    z = np.arange(0., 22000., 50.)
    u = np.zeros(len(z))
    v = np.zeros(len(z))

    # parameters
    zdep1 = self.getOption('z_curve,top')
    umax = self.getOption('u_max')
    vmax = self.getOption('v_max')
    cx = self.getOption('u_adjust')
    cy = self.getOption('v_adjust')

    for k in range(len(z)):
      if (z[k] < zdep1): # curvature section
        a = np.pi-(z[k]/zdep1)*np.pi
        u[k] = umax*np.cos(a)
        v[k] = vmax*np.sin(a)
      else: # constant section
        u[k] = umax
        v[k] = 0

    u[:] = u[:] - cx
    v[:] = v[:] - cy

    #emit
    return (z,u,v)
    
