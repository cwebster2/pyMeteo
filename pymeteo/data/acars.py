# Acars data access
# Thanks to Rich Mamrosh @ NOAA for pointing the availability of this data out to me
"""ACARS (Aircraft Communications Addressing and Reporting System) data access.

Downloads, decompresses, and parses ACARS aircraft sounding observations from
NOAA's MADIS (Meteorological Assimilation Data Ingest System) public archive.
Each sounding profile is split into ascending/descending segments and returned
with thermodynamic and kinematic variables ready for plotting on a Skew-T.
"""

import re
import netCDF4
import gzip
import io
import tempfile
import os
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from datetime import datetime
from pymeteo import dynamics, thermo
from urllib import request

data_url = "https://madis-data.ncep.noaa.gov/madisPublic1/data/point/acars/netcdf/"
airport_ids = {}


def getAvailableDatasets():
    """Fetch the list of available ACARS dataset files from MADIS.

    Scrapes the MADIS public data directory for gzip-compressed NetCDF
    filenames matching the ``YYYYMMDD_HHMM.gz`` pattern.

    :returns: List of filename strings, or *None* on failure.
    """
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
    except (OSError, ValueError) as e:
        print(f"Error fetching available datasets: {e}")


def getDataSet(dataset_name):
    """Download and decompress a single ACARS dataset from MADIS.

    :parameter dataset_name: Filename of the gzipped dataset (e.g.
        ``'20230115_1200.gz'``).
    :returns: A file-like :class:`~gzip.GzipFile` object containing the
        uncompressed NetCDF data, or *None* on failure.
    """
    req = request.Request(data_url + dataset_name)
    try:
        print("[+] Fetching dataset {0}".format(dataset_name))
        with request.urlopen(req) as f:
            compressedData = io.BytesIO(f.read())  # gzipped data
            print("[-] Decompressing response data")
            data = gzip.GzipFile(fileobj=compressedData)

            return data
    except (OSError, ValueError) as e:
        print(f"Error fetching dataset {dataset_name}: {e}")


def getAirportByCode(airport_id):
    """Look up an airport name by its numeric MADIS identifier.

    On first call the ``airport_info.dat`` file is loaded and cached in
    the module-level :data:`airport_ids` dictionary.

    :parameter airport_id: Integer airport identifier from the ACARS data.
    :returns: Airport code string (e.g. ``'KORD'``).
    :raises KeyError: If *airport_id* is not found in the data file.
    """
    print("[+] Looking up airport id '{0}'".format(airport_id))
    datfile = os.path.join(os.path.dirname(__file__), "airport_info.dat")
    if not bool(airport_ids):
        with open(datfile, "r") as f:
            for _, line in enumerate(f):
                fields = line.strip().split()
                airport_ids[int(fields[0])] = fields[1]
    return airport_ids[airport_id]


def processDataSet(data):
    """Parse a decompressed ACARS NetCDF dataset into sounding profiles.

    Reads altitude, temperature, mixing ratio, wind, and location variables
    from the NetCDF file.  Computes derived quantities (pressure, potential
    temperature, u/v wind components) and splits the data into individual
    sounding profiles based on timestamp changes.

    :parameter data: A file-like object (from :func:`getDataSet`) containing
        uncompressed NetCDF data.
    :returns: List of dicts, each representing one sounding profile with keys:
        ``'i'``, ``'n'``, ``'z'``, ``'p'``, ``'th'``, ``'qv'``, ``'u'``,
        ``'v'``, ``'lat'``, ``'lon'``, ``'airport'``, ``'time'``, ``'flag'``.
    """
    print("[+] Writing data into temporary file")
    tdata = tempfile.NamedTemporaryFile()
    tdata.write(data.read())
    print("[-] Data written to {0}".format(tdata.name))
    print("[+] Opening data as NetCDF")
    d = data.read()
    with netCDF4.Dataset(tdata.name, mode="r") as nc:
        print("[-] Dataset open with")

        _z = nc["altitude"][:]
        _T = nc["temperature"][:]
        _qv = nc["waterVaporMR"][:]
        windSpeed = nc["windSpeed"][:]
        windDir = nc["windDir"][:]
        _lon = nc["longitude"][:]
        _lat = nc["latitude"][:]
        flag = nc["sounding_flag"][:]
        _airport = nc["sounding_airport_id"][:]
        time = nc["soundingSecs"][:]

        print("[-] {0} Records".format(len(_z)))
        # conversions
        _p = thermo.p_from_pressure_altitude(_z, _T)
        _u, _v = dynamics.wind_deg_to_uv(windDir, windSpeed)
        _th = thermo.theta(_T, _p)

        # split the arrays when the flag changes sign
        splits = np.where(np.diff(time))[0] + 1

        _z = np.split(_z, splits)
        _p = np.split(_p, splits)
        _th = np.split(_th, splits)
        _qv = np.split(_qv, splits)
        _u = np.split(_u, splits)
        _v = np.split(_v, splits)
        _lat = np.split(_lat, splits)
        _lon = np.split(_lon, splits)
        _airport = np.split(_airport, splits)
        time = np.split(time, splits)
        flag = np.split(flag, splits)

        print("[-] Found {0} profiles".format(len(_z)))

        # re-shape data
        outputData = []
        for i in range(len(_z)):
            ts = time[i].compressed()
            if len(ts) == 0:
                # profiles without timestamps invalid?
                continue

            profileDir = flag[i][0]
            if profileDir == 0:
                continue

            z = _z[i].filled()
            p = _p[i].filled()
            th = _th[i].filled()
            qv = _qv[i].filled()
            u = _u[i].filled()
            v = _v[i].filled()
            lat = _lat[i].filled()
            lon = _lon[i].filled()
            try:
                airport = getAirportByCode(_airport[i][0])
            except KeyError:
                print("Airport not found ", _airport[i][0])
                airport = "ERROR"

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
                "airport": airport,
                "time": datetime.utcfromtimestamp(ts.mean()).strftime("%H%MZ"),
                "flag": profileDir,
            }
            outputData.append(profileData)

        return outputData
