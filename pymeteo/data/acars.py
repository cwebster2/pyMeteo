# Acars data access
# Thanks to Rich Mamrosh @ NOAA for pointing the availability of this data out to me

import re
import netCDF4
import gzip
import io
import tempfile
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from pymeteo import dynamics, thermo
import pymeteo.cm1.PlotWidget as PW
import pymeteo.cm1.StatWidget as SW
try:
    # For Python 3.0 and later
    from urllib import request
except ImportError:
    # Fall back to Python 2's urllib2
    import urllib2 as request

data_url = "https://madis-data.ncep.noaa.gov/madisPublic1/data/point/acars/netcdf/"

def getAvailableDatasets():
    # crawl links at https://madis-data.ncep.noaa.gov/madisPublic1/data/point/acars/netcdf/
    req = request.Request(data_url)
    linkMatcher = re.compile(r"\"([0-9_]+\.gz)\"")
    try: 
        print("[+] Fetching list of resources available")
        with request.urlopen(req) as f:
            print("[-] Parsing list")
            data = str(f.read())
            sets = linkMatcher.findall(data)
            return sets
    except:
        print("error")

def getDataSet(set):
    req = request.Request(data_url + set)
    try:
        print("[+] Fetching dataset {0}".format(set))
        with request.urlopen(req) as f:
            compressedData = io.BytesIO(f.read()) # gzipped data
            print("[-] Decompressing response data")
            data = gzip.GzipFile(fileobj=compressedData)
            
            return data
    except:
        print("error")

def processDataSet(data):
    print("[+] Writing data into temporary file")
    tdata = tempfile.NamedTemporaryFile()
    tdata.write(data.read())
    print("[-] Data written to {0}".format(tdata.name))
    print("[+] Opening data as NetCDF")
    d = data.read()
    with netCDF4.Dataset(tdata.name, mode='r') as nc:
        print("[-] Dataset open with")

        z = nc["altitude"][:]
        T = nc["temperature"][:]
        qv = nc["waterVaporMR"][:]
        windSpeed = nc["windSpeed"][:]
        windDir = nc["windDir"][:]
        lon = nc["longitude"][:]
        lat = nc["latitude"][:]
        flag = nc["sounding_flag"][:]
        airport = nc["sounding_airport_id"][:]
        time = nc["soundingSecs"][:]

        rS = nc["destAirport"]
        print ("[-] {0} Records".format(len(z)))

        #conversions

        p = thermo.p_from_pressure_altitude(z, T)
        u,v = dynamics.wind_deg_to_uv(windDir, windSpeed)
        th = thermo.theta(T, p)

        # split the arrays when the flag changes sign

        splits = np.where(np.diff(flag))[0]+1
        print(splits)

        z = np.split(z, splits)
        p = np.split(p, splits)
        th = np.split(th, splits)
        qv = np.split(qv, splits)
        u = np.split(u, splits)
        v = np.split(v, splits)
        lat = np.split(lat, splits)
        lon = np.split(lon, splits)
        airport = np.split(airport, splits)
        time = np.split(time, splits)
        flag = np.split(flag, splits)

        print("[-] Found {0} profiles".format(len(flag)))
        return (z, p, th, qv, u, v, lat, lon, airport, time, flag)        

class MainWindow(QtWidgets.QWidget):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()

    def initUI(self):
        self.splitter_2 = QtWidgets.QSplitter(QtCore.Qt.Horizontal)

        self.LVert_widget = QtWidgets.QWidget(self.splitter_2)
        self.verticalLayout = QtWidgets.QVBoxLayout(self.LVert_widget)
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

        self.label2 = QtWidgets.QLabel('Select ACARS data')
        self.verticalLayout.addWidget(self.label2)
        self.acarsFiles = QtWidgets.QLabel('Dropdown goes here')
        self.verticalLayout.addWidget(self.acarsFiles)

        self.line_2 = QtWidgets.QFrame(self.LVert_widget)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.verticalLayout.addWidget(self.line_2)

        self.label3 = QtWidgets.QLabel('Select Profile')
        self.verticalLayout.addWidget(self.label3)
        self.profile = QtWidgets.QLabel('Dropdown goes here')
        self.verticalLayout.addWidget(self.profile)

        self.drawButton = QtWidgets.QPushButton("Update Plot", self.LVert_widget)
        self.quitButton = QtWidgets.QPushButton("Quit", self.LVert_widget)
        self.verticalLayout.addWidget(self.drawButton)
        self.verticalLayout.addWidget(self.quitButton)
        self.splitter = QtWidgets.QSplitter(self.splitter_2)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.Hodograph = PW.PlotWidget(self.splitter)
        self.SoundingInfo = SW.StatWidget(self.splitter)
        self.Sounding = PW.PlotWidget(self.splitter_2)

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
        self.setWindowTitle('ACARS Profile Plotting Tool')

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

       datasets = getAvailableDatasets()
       data = getDataSet(datasets[15])
       (z, p, th, qv, u, v, lat, lon, airport, time, flag) = processDataSet(data) 
       self.Sounding.plot_sounding(z[0],th[0],p[0],qv[0],u[0],v[0])
       self.Hodograph.plot_hodograph(z[0],u[0],v[0])

    @pyqtSlot()
    def enableButton(self):
       self.drawButton.setEnabled(True)    

if __name__ == '__main__':
    datasets = getAvailableDatasets()
    data = getDataSet(datasets[15])
    (z, p, th, qv, u, v, lat, lon, airport, time, flag) = processDataSet(data) 

