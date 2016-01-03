# uwyo.py

import numpy as np
import urllib.request as request
import io

def fetch_from_web(year, month, day, hour, station):
    # http://weather.uwyo.edu/cgi-bin/sounding
    # region=naconf
    # TYPE=TEXT-LIST
    # YEAR=2016
    # MONTH=01
    # FROM=0212
    # TO=0212
    # STNM=72251

    print(year, month, day, hour, station)
    base_url = "http://weather.uwyo.edu/cgi-bin/sounding"
    url = "{0}?TYPE=TEXT%3ALIST&YEAR={1}&MONTH={2:02d}&FROM={3:02d}{4:02d}&TO={3:02d}{4:02d}&STNM={5}".format(
           base_url, year, month, day, hour, station)
    print(url)
    data= io.StringIO()
    bulkdata=[]
    urlreq = request.Request(url, method='GET')
    with request.urlopen(urlreq) as f:
        bulkdata = str(f.read()).split(r'\n')[0:-1]
        for i in range(len(bulkdata)):
            bulkdata[i] = bulkdata[i][0:]
    print(f.status)
    print(f.reason)
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
            data.write(str(bulkdata[i]) + '\n')
        elif parse_state == 'finished':
            pass

    p, z, qv, wind_dir, wind_speed, th = np.genfromtxt(io.BytesIO(data.getvalue().encode()), unpack=True,
                                                       skip_header=5,
                                                       delimiter=7,
                                                       usecols=(0,1,5,6,7,8))

    #print(title,p,z,qv,wind_dir, wind_speed, th)
