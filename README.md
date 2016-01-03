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

cm1_geninit -- visualize analytical skewt and hodograph and write sounding file suitable for CM1 or WRF  

[TOC]

#CM1 / WRF intput sounding generation

The utility:

```
$ cm1_geninit
```

allows the previwing of analytic soundings and hodographs and writing to files suitable for initializing
the CM1 or WRF weather models.  The code is extensible and adding your own analytic profiles should be
straightforward.  Future documentation will detail this and I welcome pull requests!

#Skew-T / Ln-P diagrams

Skew-T / Ln-P plots are simply 2-D plots with a skewed temperature axis and a logarithmic pressure axes (y).  Data plotted on this style of plot include temperature, dew point, and various other derived values such as the temperature of a parcel lifted from the surface.  The plots produced by these scripts also have a hodograph plotted and a data block with convective storm parameters.

## Example plot

[![skewt](https://wxster.com/static/media/skewt/skewt.png)](https://wxster.com/static/media/skewt/skewt.png)

##Plotting from data

In all of the plotting methods, the plot output type is determined by the extension of the output file provided.  You can write any type of file that the matplotlib backend can write.

###Tabular

If you have tabular data suitable for WRF or CM1 model initialization, you can plot a skewt of this data with:

```
$ skewt tabular -f sounding.dat skewt.pdf
```

The format of the sounding data file is

> 1 line header that contains:  surface pressure (mb)    surface theta (K)    surface qv (g/kg)
>
> following lines are:  height (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)

See the file `testdata/sounding_wrfinit.dat` for an example of this file format.

###CM1

For CM1 output in native GrADS format, you can plot a skewt from model output with

```
$ skewt cm1 -p . -d cm1out -x 0 -y 0 skewt-cm1.pdf
```
In this case, `-p` is the path to the dataset, `-d` is the CM1 `output_basename`, `-x` and `-y` are the location of the plot in km and `-o` is the file to output.  The current version (v0.4) of this script only works for CM1 datasets that are output with one file per timestep and will plot whatever timestep the file contains. 

###CM1 / HDF5

For CM1 output in HDF5 format, you can plot a skewt from model output with:

```
$ skewt cm1hdf -f model-data.h5 -x 0 -y 0 output.pdf
```

This currently requires that HDF5 output be un-tiled (though it might work for tiled files, I have not checked).  The options are as in the CM1 GrADS version except that there is no `-p` and `-f` references a specific HDF5 file.  This also assumes one file per timestep and will plot whatever timestep the file contains.

###WRF

For WRF output in NetCDF format, you can plot a skewt from model output with:

```
$ skewt wrf -d wrfou.nc --lat 30 --lon -80 -t 0 skewt.pdf
```

In this case, `-f` references a WRF output file, `--lat` and `--lon` reference a location within the WRF domain, `-t` reference a timestep within the WRF output and `-o` specifies an output files.

####University of Wyoming sounding data

From a file:

```
$ skewt uwyo -f uwyo-data.dat skewt.pdf
```

From the website for the most recent sounding from a station: 

```
$ skewt uwyoweb --station 72251 skewt.pdf
```

##From Numpy arrays

```
import numpy as np
import pymeteo.skewt as skewt

# prepare 1D arrays height (z), pressure (p), potential temperature (th), 
# water vapor mixing ratio (qv), winds (u and v) all of the same length.

skewt.plot(None, z, th, p, qv, u, v, 'output,pdf')
```

You can also choose to plot just the sounding, just the hodograph or plot each on on axes that you define.  For details see the implementation of `pymeteo.skewt.plot()` and the [pyMeteo documentation][1].


  [1]: http://pythonhosted.org/pymeteo/