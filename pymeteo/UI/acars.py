from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from pymeteo.data import acars
import pymeteo.cm1.PlotWidget as PW
import pymeteo.cm1.StatWidget as SW

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

        self.acarsFiles = QtWidgets.QComboBox(self)
        self.verticalLayout.addWidget(self.acarsFiles)

        self.line_2 = QtWidgets.QFrame(self.LVert_widget)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.verticalLayout.addWidget(self.line_2)

        self.label3 = QtWidgets.QLabel('Select Profile')
        self.verticalLayout.addWidget(self.label3)
        self.profiles = QtWidgets.QComboBox(self)
        self.verticalLayout.addWidget(self.profiles)

        self.verticalLayout.addStretch(1)

        self.drawButton = QtWidgets.QPushButton("Update Plot", self.LVert_widget)
        self.quitButton = QtWidgets.QPushButton("Quit", self.LVert_widget)
        self.verticalLayout.addWidget(self.drawButton)
        self.verticalLayout.addWidget(self.quitButton)
        self.splitter = QtWidgets.QSplitter(self.splitter_2)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.Hodograph = PW.PlotWidget(self.splitter)
        self.SoundingInfo = SW.StatWidget(self.splitter)
        self.Sounding = PW.PlotWidget(self.splitter_2)

        datasets = acars.getAvailableDatasets()
        self.acarsFiles.addItems(datasets)
        self.acarsFiles.activated[str].connect(self.processDataset)

        self.profiles.activated.connect(self.displayProfile)

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

    def profileDesc(self,profile):
        profileDir = "UP" if profile["flag"] == 1 else "DOWN"
        return "{0}: {1}({2}) {3}".format(profile["i"], profile["airport"], profileDir, profile["time"])

    def processDataset(self, acarsFile):
        print("[+] Getting Data")
        self.profiles.clear()
        self.data = acars.processDataSet(acars.getDataSet(acarsFile)) 
        for profile in self.data:
            desc = self.profileDesc(profile)
            print("[-] Adding profile {0}".format(desc))
            self.profiles.addItem(desc)

    def displayProfile(self, index):
        print("[+] Displaying profile {0}".format(index))
        profile = self.data[index]
        desc = self.profileDesc(profile)
        print("[-] {0}".format(desc))
        self.Sounding.plot_sounding(
            profile["z"],
            profile["th"],
            profile["p"],
            profile["qv"],
            profile["u"],
            profile["v"]
        )
        self.Hodograph.plot_hodograph(
            profile["z"],
            profile["u"],
            profile["v"]
        )


    exportSounding = pyqtSignal(str)

    def export(self):
       self.exportSounding.emit(self.lineEdit.text())

    def update_plot(self):
       self.drawButton.setEnabled(False)


    @pyqtSlot()
    def enableButton(self):
       self.drawButton.setEnabled(True)    
