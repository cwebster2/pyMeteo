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
from datetime import datetime
from pymeteo import dynamics, thermo
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

        _z = nc["altitude"][:]
        _T = nc["temperature"][:]
        _qv = nc["waterVaporMR"][:]
        windSpeed = nc["windSpeed"][:]
        windDir = nc["windDir"][:]
        _lon = nc["longitude"][:]
        _lat = nc["latitude"][:]
        flag = nc["sounding_flag"][:]
        airport = nc["sounding_airport_id"][:]
        time = nc["soundingSecs"][:]

        rS = nc["destAirport"] #flight, tailnumber
        print ("[-] {0} Records".format(len(_z)))
        print(rS)
        #conversions
        _p = thermo.p_from_pressure_altitude(_z, _T)
        _u, _v = dynamics.wind_deg_to_uv(windDir, windSpeed)
        _th = thermo.theta(_T, _p)

        # split the arrays when the flag changes sign
        splits = np.where(np.diff(flag))[0]+1

        _z = np.split(_z, splits)
        _p = np.split(_p, splits)
        _th = np.split(_th, splits)
        _qv = np.split(_qv, splits)
        _u = np.split(_u, splits)
        _v = np.split(_v, splits)
        _lat = np.split(_lat, splits)
        _lon = np.split(_lon, splits)
        airport = np.split(airport, splits)
        time = np.split(time, splits)
        flag = np.split(flag, splits)

        for i in range(len(_z[0])):
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(_z[0][i], _p[0][i], _th[0][i], _qv[0][i],
            _u[0][i], _v[0][i], flag[0][i]))

        print("[-] Found {0} profiles".format(len(_z)))
        #re-shape data
        outputData = []
        for i in range(len(_z)):
            ts = time[i].compressed()
            if len(ts) == 0:
                # profiles without timestamps invalid?
                continue 

            profileDir = flag[i][0]
            z = _z[i].filled()
            p = _p[i].filled()
            th = _th[i].filled()
            qv = _qv[i].filled()
            u = _u[i].filled()
            v = _v[i].filled()
            lat = _lat[i].filled()
            lon = _lon[i].filled()
            profileData = {
                "i": i,
                "n": len(z),
                "z": z if profileDir > 0 else z[::-1],
                "p": p if profileDir > 0 else p[::-1],
                "th": th if profileDir > 0 else th[::-1],
                "qv": qv if profileDir > 0 else qv[::-1],
                "u": u if profileDir > 0 else u[::-1],
                "v": v if profileDir > 0 else v[::-1],
                "lat": lat if profileDir > 0 else lat[::-1],
                "lon": lon if profileDir > 0 else lon[::-1],
                "airport": airport[i][0],
                "time": datetime.utcfromtimestamp(ts.mean()).strftime("%H%MZ"),
                "flag": profileDir
            }
            outputData.append(profileData)

        return outputData
