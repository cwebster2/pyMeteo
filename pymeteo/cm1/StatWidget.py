from __future__ import print_function
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSlot, pyqtSignal
import numpy as np

from pymeteo import skewt
from pymeteo import dynamics

# widgets -> skewt wind sounding etc

class StatWidget(QtWidgets.QFrame):

  statsdone = pyqtSignal()

  def __init__(self, parent=None):
    super(StatWidget, self).__init__(parent)
    self.initUI()
    self.z = 0
    self.u = 0
    self.v = 0
    self.th = 0
    self.p = 0
    self.qv = 0
    self.SStats = 0
    self.HStats = 0

  def initUI(self):
    self.statsBox = QtWidgets.QTextBrowser(self)
    #self.updateBtn = QtGui.QPushButton("Update", self)
    #self.updateBtn.clicked.connect(self.calcStats)

    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(self.statsBox)
    #layout.addWidget(self.updateBtn)
    self.setLayout(layout)

  def calcStats(self):
    self.statsBox.clear()
    if (self.SStats == 0):
      self.statsBox.append("No Sounding Plotted")  
    else:
      pcl, mupcl, mlpcl = skewt.calc_sounding_stats(self.z, self.th, self.p, self.qv)
      self.statsBox.append("Var\t\tSFC\tML\tMU")
      self.statsBox.append("CAPE\tJ/kg\t{0:.0f}\t{1:.0f}\t{2:.0f}".format(pcl['cape'],mlpcl['cape'],mupcl['cape']))
      self.statsBox.append("CIN\tJ/kg\t{0:.1f}\t{1:.1f}\t{2:.1f}".format(pcl['cin'],mlpcl['cin'],mupcl['cin']))
      self.statsBox.append("LI_max\t\t{0:.1f}\t{1:.1f}\t{2:.1f}".format(pcl['max_li'],mlpcl['max_li'],mupcl['max_li']))
      self.statsBox.append("LI_500\t\t{0:.1f}\t{1:.1f}\t{2:.1f}".format(pcl['li500'],mlpcl['li500'],mupcl['li500']))
      self.statsBox.append("LI_300\t\t{0:.1f}\t{1:.1f}\t{2:.1f}".format(pcl['li300'],mlpcl['li300'],mupcl['li300']))
      self.statsBox.append("TOPS\tm\t{0:.0f}\t{1:.0f}\t{2:.0f}".format(pcl['ztops'],mlpcl['ztops'],mupcl['ztops']))

      self.statsBox.append("LCL\tmb (m)\t{0:.0f} ({1:.0f})\t{2:.0f} ({3:.0f})\t{4:.0f} ({5:.0f})".format(pcl['lclprs']/100.,pcl['lcl'],mlpcl['lclprs']/100.,mlpcl['lcl'],mupcl['lclprs']/100.,mupcl['lcl']))
      self.statsBox.append("LCL\tmb (m)\t{0:.0f} ({1:.0f})\t{2:.0f} ({3:.0f})\t{4:.0f} ({5:.0f})".format(pcl['lfcprs']/100.,pcl['lfc'],mlpcl['lfcprs']/100.,mlpcl['lfc'],mupcl['lfcprs']/100.,mupcl['lfc']))

      self.statsBox.append("PRS\tmb\t{0:.0f}\t{1:.0f}\t{2:.0f}".format(pcl['prs']/100.,mlpcl['prs']/100.,mupcl['prs']/100.))
      self.statsBox.append("")
      shear = skewt.calc_hodograph_stats(self.z, self.u, self.v)
      cth,cr = dynamics.uv_to_deg(shear['bunkers'][0],shear['bunkers'][1])
      self.statsBox.append("Storm motion (left mover): {0:.0f} deg {1:5.2f} m/s".format(cth,cr))
      self.statsBox.append("Storm motion (left mover): u={0:5.2f} v={1:5.2f} m/s".format(shear['bunkers'][0],shear['bunkers'][1]))
      self.statsBox.append("")
      self.statsBox.append("Layer\tSRH\tERH")
      self.statsBox.append("0-1 km\t{0:.0f}\t{1:.0f}".format(shear['srh01'],shear['erh01']))
      self.statsBox.append("0-3 km\t{0:.0f}\t{1:.0f}".format(shear['srh03'],shear['erh03']))
      self.statsBox.append("")
      self.statsBox.append("Layer\tShear Vector")
      self.statsBox.append("0-1 km\t{0:.0f} deg, {1:5.2f} m/s".format(*shear['s01']))
      self.statsBox.append("0-3 km\t{0:.0f} deg, {1:5.2f} m/s".format(*shear['s03']))
      self.statsBox.append("0-6 km\t{0:.0f} deg, {1:5.2f} m/s".format(*shear['s06']))

      self.statsdone.emit()

  @pyqtSlot(np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray)
  def SStat_values(self, z,th,p,qv,u,v):
    self.z = z
    self.th = th
    self.p = p
    self.qv = qv
    self.u = u
    self.v = v
    self.SStats = 1
    self.calcStats()

  @pyqtSlot(str)
  def output_sounding(self, filename):
    print("would output to file {0}".format(filename))
    f = open(filename, 'w')
    for k in range(len(self.z)):
      if (k==0):
        print("{0:8.2f}\t{1:8.6f}\t{2:8.6f}".format(self.p[0]/100., self.th[0], self.qv[0]*1000.), file=f)
      else:
        print("{0:8.2f}\t{1:8.6f}\t{2:8.6f}\t{3:8.6f}\t{4:8.6f}".format(self.z[k], self.th[k], self.qv[k]*1000., self.u[k], self.v[k]), file=f)
    f.close() 
