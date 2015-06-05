from PyQt5 import QtWidgets, QtGui, QtCore

class OptionsWidget(QtWidgets.QFrame):

  def __init__(self):
    super(OptionsWidget,self).__init__()
 
  def initUI(self, name, v):
    self.qlbl = QtWidgets.QLabel(name)
    self.grid = QtWidgets.QGridLayout()
    self.grid.addWidget(self.qlbl,1,0,1,3)

    nv = len(v)
    self.variables = {}

    for var in range(nv):
      self.variables[v[var][0]] = self.addOption(var, *v[var])

    self.grid.setSpacing(5)
    self.setLayout(self.grid)

  def addOption(self, pos, name, defval, units):
    label = QtWidgets.QLabel(name)
    var   = QtWidgets.QLineEdit(defval, self)
    label2= QtWidgets.QLabel(units)

    self.grid.addWidget(label,pos+2,0)
    self.grid.addWidget(var,pos+2,1)
    self.grid.addWidget(label2,pos+2,2)

    return var
    
  def getOption(self, name):
    return float(self.variables[name].text())

