"""
.. module: pymeteo.uwyo
   :platform: Unix, Windows
   :synopsis: University of Wyoming sounding data helper functions

.. moduleauthor:: Casey Webster <casey.webster@gmail.com>

This module provides functions for working with University of Wyoming sounding
data.  This includes functions for getting data from files and from the University
of Wyoming website at http://weather.uwyo.edu/ and a function for preparing that data
for SkewT plotting.

Getting Data
++++++++++++

* :py:func:`fetch_from_file` -- loads sounding data from a file 
* :py:func:`fetch_from_web` -- loads sounding data from website

Preparing Data
++++++++++++++

* :py:func:`transform_and_check_data` -- quality checks and converts units of data

Module Reference
++++++++++++++++
"""

import numpy as np
import io
import datetime
import pymeteo.dynamics as dynamics
try:
    # For Python 3.0 and later
    from urllib import request
except ImportError:
    # Fall back to Python 2's urllib2
    import urllib2 as request

def fetch_from_file(filename):
    """Load Uinversity of Wyoming sound data from a file

    :param filename: The filename containing sounding data
    :type filename: str
    :returns: pressure, heights, mixing ration, wind direction, wind speed and potential temperature

    This loads sounding data from a file and returns the data columns needed
    for SkewT plotting.  The data format is the same as copy/pasting data from
    the uwyo website.
    """
    with open(filename, 'r') as f:
        title = f.readline()
        skiprows = 7
        if "Obs" not in title:
            title = filenme
            skiprows = 5

    p, z, qv, wind_dir, wind_speed, th = np.genfromtxt(filename, unpack=True, skip_header=skiprows,
                                                           delimiter=(7,7,7,7,7,7,7,7,7,7,7,7),
                                                           usecols=(0,1,5,6,7,8))
    return (title, p, z, qv, wind_dir, wind_speed, th)

def fetch_from_web(date, station):
    """Load Uinversity of Wyoming sound data from the uwyo website

    :param date: Time and date of requested sounding data
    :type date: datetime
    :param station: The station ID for sounding data
    :type station: int
    :returns: pressure, heights, mixing ration, wind direction, wind speed and potential temperature

    This loads sounding data from the web and returns the data columns needed
    for SkewT plotting.  The date is rounded to the preceding 12-hourly observation at 12Z or 00Z.
    """
    # http://weather.uwyo.edu/cgi-bin/sounding
    # region=naconf
    # TYPE=TEXT-LIST
    # YEAR=2016
    # MONTH=01
    # FROM=0212
    # TO=0212
    # STNM=72251

    year = date.year
    month = date.month
    day = date.day
    hour = date.hour
    if hour < 12:
        hour = 00
    else:
        hour = 12
    hour = 00
    print(year, month, day, hour, station)
    base_url = "http://weather.uwyo.edu/cgi-bin/sounding"
    url = "{0}?TYPE=TEXT%3ALIST&YEAR={1}&MONTH={2:02d}&FROM={3:02d}{4:02d}&TO={3:02d}{4:02d}&STNM={5}".format(
           base_url, year, month, day, hour, station)
    #print(url)
    data=[]
    bulkdata=[]
    urlreq = request.Request(url) #, method='GET')

    try:
        with request.urlopen(urlreq) as f:
            bulkdata = str(f.read()).split(r'\n')[0:-1]
            for i in range(len(bulkdata)):
                bulkdata[i] = bulkdata[i][0:]

        if f.status != 200:
            print("Error fetching data");
            return

    except AttributeError:
        print("python2")
        f = request.urlopen(urlreq)
        bulkdata = f.read().split('\n')
        for i in range(len(bulkdata)):
            bulkdata[i] = bulkdata[i][0:]

    parse_state = 'start'
    for i in range(len(bulkdata)):
        if parse_state == 'start':
            if '<H2>' not in bulkdata[i]:
                continue
            title = (bulkdata[i][4:-5])
            parse_state = 'data'
        elif parse_state == 'data':
            if '</PRE>' in bulkdata[i]:
                parse_state = 'finished'
                continue
            if '<PRE>' in bulkdata[i]:
                continue
            data.append(bulkdata[i])
        elif parse_state == 'finished':
            pass
    data = "\n".join(data)

    if (data == ""):
        print("ERROR: No sounding data found in data returned from the server");
        quit();
        
    p, z, qv, wind_dir, wind_speed, th = np.genfromtxt(io.BytesIO(data.encode()), unpack=True,
                                                       skip_header=5,
                                                       delimiter=7,
                                                       usecols=(0,1,5,6,7,8))

    return (title,p,z,qv,wind_dir, wind_speed, th)


def transform_and_check_data(p, z, qv, wind_dir, wind_speed, th):
    """Quality check and convert units of Uinversity of Wyoming sound data

    :param p: pressure
    :param z: heights
    :param qv: mixing ration
    :param wind_dir: wind direction
    :param wind_speed: wind speed
    :param th: potential temperature
    :returns: pressure, heights, mixing ration, u wind, v wind and potential temperature

    This ingests data loaded from file owr web and culls invalid rows (NaN in certain colums)
    and quality checks the data.  Units are converted from the native uwyo data to the
    module standard SI units.
    """
    # clean up NaNs and mark invalid rows for deletion
    nk = len(z)
    for k in np.arange(nk):
        delete_rows = []
        if np.isnan(p[k]):
            delete_rows.append(k)
        if np.isnan(qv[k]):
            qv[k] = 0
        if np.isnan(wind_speed[k]):
            wind_speed[k] = wind_speed[k-1]
        if np.isnan(wind_dir[k]):
            wind_dir[k] = wind_dir[k-1]

    # delete invalid rows
    p = np.delete(p,delete_rows)
    z = np.delete(z,delete_rows)
    qv = np.delete(qv,delete_rows)
    wind_dir = np.delete(wind_dir, delete_rows)
    wind_speed = np.delete(wind_speed, delete_rows)
    th = np.delete(th, delete_rows)

    # convert ingested units to our package standard units
    nk = len(z)
    p = p * 100. # hPa to Pa
    qv = qv / 1000. # g/kg to kg/kg
    wind_speed = wind_speed * 0.51444  # kts to m/s

    # convert wind direction,speed to u,v components
    u = np.empty(nk, np.float32)
    v = np.empty(nk, np.float32)
    for k in np.arange(nk):
        u[k], v[k] = dynamics.wind_deg_to_uv(wind_dir[k], wind_speed[k])

    p[np.isnan(p)] = 0
    qv[np.isnan(qv)] = 0
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    #reutrn data
    return (p, z, qv, u, v, th)
