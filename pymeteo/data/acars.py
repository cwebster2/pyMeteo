# Acars data access
# Thanks to Rich Mamrosh @ NOAA for pointing the availability of this data out to me

import re
import netCDF4
import gzip
import io
import tempfile
import numpy as np
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

        # split the arrays when the flag changes sign

        splits = np.where(np.diff(flag))[0]+1
        print(splits)

        z = np.split(z, splits)
        p = np.split(p, splits)
        T = np.split(T, splits)
        qv = np.split(qv, splits)
        u = np.split(u, splits)
        v = np.split(v, splits)
        lat = np.split(lat, splits)
        lon = np.split(lon, splits)
        airport = np.split(airport, splits)
        time = np.split(time, splits)
        flag = np.split(flag, splits)

        print("[-] Found {0} profiles".format(len(flag)))
        return (z, p, T, qv, u, v, lat, lon, airport, time, flag)        
            
if __name__ == '__main__':
    datasets = getAvailableDatasets()
    data = getDataSet(datasets[15])
    processDataSet(data)

