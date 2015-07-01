# pyMeteo
General meteorological routines, skew-T/log-p plotting and working with CM1 model data.

Online documentation available at https://pythonhosted.org/pymeteo/

## Getting pyMeteo

The easiest way to get pymeteo is via pip

    pip install pymeteo

pyMeteo is developed with python 3.4 and should work with any version of python 3 but may not work
with python 2.  Open an issue or send me a pull request if you want to make python 2 work but make
sure any submissions do not regress against python 3.

You can obtain the most recent development version of pymeteo by cloing this repository.

    git clone https://github.com/cwebster2/pyMeteo.git
    cd pyMeteo
    python3 setup.py install

## Inlcuded scripts

skewt -- plots skewt from CM1 output (Grads format) or from CM1 or WRF input soundings  
skewt-hdf -- plots skewt from CM1 output (HDF5 format, processed by ingest_cm1)  
cm1_geninit -- visualize analytical skewt and hodograph and write sounding file suitable for CM1 or WRF  

