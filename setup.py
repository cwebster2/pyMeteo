"""Meteorology utilities, SkewT plotting and more

This package contains some general meteorological routines 
and some specific routines to work with CM1.  The general
routines include generation of Skew-T/Log-P plotting.  The CM1
routines include generation of input soundings for CM1 and
classes for reading CM1 output files in grads and HDF5 format.
"""

classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Topic :: Scientific/Engineering :: Atmospheric Science
Development Status :: 3 - Alpha
Topic :: Software Development :: Libraries :: Python Modules
Programming Language :: Python :: 3
"""

from distutils.core import setup
#from setuptools import setup

def readme():
  with open('README.md') as f:
    return f.read()

doclines = __doc__.split("\n")

setup(name='pymeteo',
      version='0.4.0',
      description=doclines[0],
      long_description="\n".join(doclines[2:]),
      url='http://github.com/cwebster2/pymeteo',
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
               'bin/skewt-blank'],
      classifiers=filter(None, classifiers.split("\n")))

