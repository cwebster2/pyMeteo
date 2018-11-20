#!/usr/bin/env python

import sys, os
import matplotlib
matplotlib.use('qt5agg')
from pyQt5 import QtWidgets, QtGui, QtCore
from pymeteo.data.acars import MainInterface

def main():
    app = QtWidgets.QApplication(sys.argv)
    w = MainInterface.MainWindow()
    sys.exist(app.exec_())

if __name__ == '__main__':
    main()