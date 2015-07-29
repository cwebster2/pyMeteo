from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as \
    FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as \
    NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np

from pymeteo.cm1.soundings import WK82
from pymeteo import skewt

# widgets -> skewt wind sounding etc


class PlotWidget(QtWidgets.QFrame):

    SStats = pyqtSignal(np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                        np.ndarray, np.ndarray)

    def __init__(self, parent=None):
        super(PlotWidget, self).__init__(parent)
        self.initUI()

    def initUI(self):
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot_sounding_axes(self):
        self.figure.clear()
        ax = self.figure.add_axes([0.05, 0.05, 0.90, 0.945])
        # ax = self.figure.add_axes([0.005,0.05,0.985,0.945])
        ax.hold(True)
        skewt.set_fontscalefactor(4)
        skewt.plot_sounding_axes(ax)
        self.canvas.draw()

    def plot_hodograph_axes(self):
        ax = self.figure.add_axes([0.005, 0.005, 0.985, 0.985])
        ax.hold(True)
        skewt.plot_hodo_axes(ax)
        self.canvas.draw()

    def plot_sounding(self, z, th, p, qv, u, v):
        self.figure.clear()
        # ax = self.figure.add_axes([0.005,0.05,0.985,0.945])
        ax = self.figure.add_axes([0.05, 0.05, 0.90, 0.945])
        ax.hold(True)
        skewt.plot_sounding_axes(ax)
        skewt.plot_sounding(ax, z, th, p, qv, u, v)
        self.canvas.draw()
        # Send data to stats widget
        self.SStats.emit(z, th, p, qv, u, v)

    def plot_hodograph(self, z, u, v):
        self.figure.clear()
        ax = self.figure.add_axes([0.005, 0.05, 0.985, 0.945])
        ax.hold(True)
        skewt.plot_hodo_axes(ax)
        skewt.plot_hodograph(ax, z, u, v)
        self.canvas.draw()
