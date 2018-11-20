# Acars data access
# Thanks to Rich Mamrosh @ NOAA for pointing the availability of this data out to me

import re
import netCDF4
import gzip
import io
import tempfile
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

        altitude = nc["altitude"]
        temperature = nc["temperature"]
        dewpoint = nc["dewpoint"]
        windSpeed = nc["windSpeed"]
        windDir = nc["windDir"]
        lon = nc["longitude"]
        lat = nc["latitude"]
        flag = nc["sounding_flag"]
        airport = nc["sounding_airport_id"]
        time = nc["soundingSecs"]

        rS = nc["destAirport"]
        print ("[-] {0} Records".format(len(altitude)))
        for i in range(100):
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(
                altitude[i], temperature[i], dewpoint[i], windSpeed[i], windDir[i],
                lat[i], lon[i], flag[i], airport[i], time[i], rS[i,:]
            ))
        
        # vars: longitude, latitude, sounding_flag, sounding_airport_id
        # z altitude (m)
        # th  temperature (K)
        # p altitude 
        # qv dewpoint (K)
        # u - winDir windSpeed m/s
        # v
        # sounding_secs
            
if __name__ == '__main__':
    datasets = getAvailableDatasets()
    data = getDataSet(datasets[15])
    processDataSet(data)

