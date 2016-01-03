# uwyo.py

import numpy as np
import urllib.request as request
import io
import datetime
import pymeteo.dynamics as dynamics


def fetch_from_file(filename):
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
    print(year, month, day, hour, station)
    base_url = "http://weather.uwyo.edu/cgi-bin/sounding"
    url = "{0}?TYPE=TEXT%3ALIST&YEAR={1}&MONTH={2:02d}&FROM={3:02d}{4:02d}&TO={3:02d}{4:02d}&STNM={5}".format(
           base_url, year, month, day, hour, station)
    #print(url)
    data=[]
    bulkdata=[]
    urlreq = request.Request(url, method='GET')
    with request.urlopen(urlreq) as f:
        bulkdata = str(f.read()).split(r'\n')[0:-1]
        for i in range(len(bulkdata)):
            bulkdata[i] = bulkdata[i][0:]
    if f.status != 200:
        print("Error fetching data");
        return
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

    p, z, qv, wind_dir, wind_speed, th = np.genfromtxt(io.BytesIO(data.encode()), unpack=True,
                                                       skip_header=5,
                                                       delimiter=7,
                                                       usecols=(0,1,5,6,7,8))

    return (title,p,z,qv,wind_dir, wind_speed, th)


def transform_and_check_data(p, z, qv, wind_dir, wind_speed, th):
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

    #reutrn data
    return (p, z, qv, u, v, th)
