"""This module provides constants used in the rest of the package

"""
# some constants

missingval = -99999999.
m2km = 0.001
km2m = 1000.
gravity = 9.81
maxparcels = 99999
L = 2.501e6    # latent heat of vaporization
Rd = 287.04         # gas constant dry air
Rv = 461.5              # gas constant water vapor
epsilon = Rd/Rv
cp = 1005.7              # what about cpd vs cpv
cpd = 1005.7             # what about cpd vs cpv
cpv = 1875.0
cpl = 4190.0
cpi = 2118.636
cv = 718.
g = 9.81
p00 = 100000.   # reference pressure
T00 = 273.15
xlv = L
xls = 2836017.0

# Derivced values

lv1 = xlv+(cpl-cpv)*T00
lv2 = cpl - cpv
ls1 = xls+(cpi-cpv)*T00
ls2 = cpi - cpv

kappa = (cp-cv)/cp
kappa_d = Rd/cp
rp00 = 1./p00
reps = Rv/Rd
eps = epsilon
rddcp = kappa_d
cpdrd = cp/Rd
cpdg = cp/g

converge = 0.0002

