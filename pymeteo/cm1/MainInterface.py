import sys
import inspect
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot

import pymeteo.cm1.soundings
import pymeteo.cm1.hodographs

from pymeteo import skewt
import pymeteo.cm1.PlotWidget as PW
import pymeteo.cm1.StatWidget as SW


class MainWindow(QtWidgets.QWidget):

   def __init__(self):
      super(MainWindow,self).__init__()
      self.initUI()

   def initUI(self):
      
      self.splitter_2 = QtWidgets.QSplitter(QtCore.Qt.Horizontal)

      self.LVert_widget = QtWidgets.QWidget(self.splitter_2)
      self.verticalLayout = QtWidgets.QVBoxLayout(self.LVert_widget)
      #self.verticalLayout.setMargin(0)
      self.horizontalLayout = QtWidgets.QHBoxLayout()
      self.label = QtWidgets.QLabel('Output file:',self.LVert_widget)
      self.horizontalLayout.addWidget(self.label)
      self.lineEdit = QtWidgets.QLineEdit(self.LVert_widget)
      self.horizontalLayout.addWidget(self.lineEdit)
      self.pushButton = QtWidgets.QPushButton("Export", self.LVert_widget)
      self.horizontalLayout.addWidget(self.pushButton)
      self.verticalLayout.addLayout(self.horizontalLayout)
      self.line = QtWidgets.QFrame(self.LVert_widget)
      self.line.setFrameShape(QtWidgets.QFrame.HLine)
      self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
      self.verticalLayout.addWidget(self.line)
      self.SoundingTabs = QtWidgets.QTabWidget(self.LVert_widget)
      self.verticalLayout.addWidget(self.SoundingTabs)
      self.line_2 = QtWidgets.QFrame(self.LVert_widget)
      self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
      self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
      self.verticalLayout.addWidget(self.line_2)
      self.WindTabs = QtWidgets.QTabWidget(self.LVert_widget)
      self.verticalLayout.addWidget(self.WindTabs)
      self.drawButton = QtWidgets.QPushButton("Update Plot", self.LVert_widget)
      self.quitButton = QtWidgets.QPushButton("Quit", self.LVert_widget)
      self.verticalLayout.addWidget(self.drawButton)
      self.verticalLayout.addWidget(self.quitButton)
      self.splitter = QtWidgets.QSplitter(self.splitter_2)
      self.splitter.setOrientation(QtCore.Qt.Vertical)
      self.Hodograph = PW.PlotWidget(self.splitter)
      self.SoundingInfo = SW.StatWidget(self.splitter)
      self.Sounding = PW.PlotWidget(self.splitter_2)

      # Load all of the sounding modules into tabs
      for name,obj in inspect.getmembers(pymeteo.cm1.soundings):
        if inspect.ismodule(obj):
          cobj = getattr(obj,name)
          if inspect.isclass(cobj):
             instance = cobj()
             self.SoundingTabs.addTab(instance, name)
      
      # Load all of the hodograph modules into tabs
      for name,obj in inspect.getmembers(pymeteo.cm1.hodographs):
        if inspect.ismodule(obj):
          cobj = getattr(obj,name)
          if inspect.isclass(cobj):
             instance = cobj()
             self.WindTabs.addTab(instance, name)

      self.Sounding.SStats.connect(self.SoundingInfo.SStat_values)

      self.pushButton.clicked.connect(self.export)
      self.exportSounding.connect(self.SoundingInfo.output_sounding)

      self.drawButton.clicked.connect(self.update_plot)
      self.quitButton.clicked.connect(self.close)

      self.hbox = QtWidgets.QHBoxLayout(self)
      self.hbox.addWidget(self.splitter_2)
      self.setLayout(self.hbox)

      self.resize(1200,768)

      #TODO: size widgets initially. This doesnt work
      self.LVert_widget.resize(250,500)
      self.Hodograph.resize(250,250)
      self.SoundingInfo.resize(250,250)
      self.Sounding.resize(250,500)


      self.center()
      self.setWindowTitle('CM1 / WRF sounding initialization tool')

      self.show()
      self.Sounding.plot_sounding_axes()
      self.Hodograph.plot_hodograph_axes()
      self.SoundingInfo.calcStats()

      self.SoundingInfo.statsdone.connect(self.enableButton)

   def center(self):
      qr = self.frameGeometry()
      cp = QtWidgets.QDesktopWidget().availableGeometry().center()
      qr.moveCenter(cp)
      self.move(qr.topLeft())

   exportSounding = pyqtSignal(str)

   def export(self):
     self.exportSounding.emit(self.lineEdit.text())

   def update_plot(self):
     self.drawButton.setEnabled(False)
     z,th,p,qv = self.SoundingTabs.currentWidget().plot()
     z2,u,v = self.WindTabs.currentWidget().plot()
     self.Sounding.plot_sounding(z,th,p,qv,u,v)
     self.Hodograph.plot_hodograph(z2,u,v)

   @pyqtSlot()
   def enableButton(self):
     self.drawButton.setEnabled(True)
