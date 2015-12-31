"""Meteorology utilities, SkewT plotting and more

This package contains some general meteorological routines 
and some specific routines to work with CM1 and WRF.  The general
routines include generation of Skew-T/Log-P plotting.  The CM1
routines include generation of input soundings for CM1 and
classes for reading CM1 output files in grads and HDF5 format.
"""

classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Topic :: Scientific/Engineering :: Atmospheric Science
Development Status :: 4 - Beta
Topic :: Software Development :: Libraries :: Python Modules
Programming Language :: Python :: 3
Programming Language :: Python :: 3.4
Programming Language :: Python :: 2
Programming Language :: Python :: 2.7
"""

import os
import subprocess
from distutils.core import setup
#from setuptools import setup

def readme():
  with open('README.md') as f:
    return f.read()

doclines = __doc__.split("\n")

version_py = os.path.join(os.path.dirname(__file__), 'VERSION')
with open(version_py, 'r') as fh:
    version_git = open(version_py).read().strip()

setup(name='pymeteo',
      version='{ver}'.format(ver=version_git),
      description=doclines[0],
      long_description="\n".join(doclines[2:]),
      url='https://wxster.com/projects/pymeteo',
      maintainer="Casey Webster",
      maintainer_email="casey.webster@gmail.com",
      author='Casey Webster',
      author_email='casey.webster@gmail.com',
      license='BSD 3-clause',
      platforms= ["any"],
      packages=['pymeteo','pymeteo.cm1',
                'pymeteo.cm1.hodographs',
                'pymeteo.cm1.soundings'],
      scripts=['bin/cm1_geninit',
               'bin/skewt',
               'bin/skewt-hdf',
               'bin/skewt-blank',
               'bin/skewt-wrf'],
      classifiers=filter(None, classifiers.split("\n")))

