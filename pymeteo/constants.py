"""This module provides constants used in the rest of the package"""
# some constants

missingval: float = -99999999.0
m2km: float = 0.001
km2m: float = 1000.0
gravity: float = 9.81
maxparcels: int = 99999
L: float = 2.501e6  # latent heat of vaporization
Rd: float = 287.04  # gas constant dry air
Rv: float = 461.5  # gas constant water vapor
epsilon: float = Rd / Rv
cp: float = 1005.7  # what about cpd vs cpv
cpd: float = 1005.7  # what about cpd vs cpv
cpv: float = 1875.0
cpl: float = 4190.0
cpi: float = 2118.636
cv: float = 718.0
g: float = 9.81
p00: float = 100000.0  # reference pressure
T00: float = 273.15
xlv: float = L
xls: float = 2836017.0

# Derived values

lv1: float = xlv + (cpl - cpv) * T00
lv2: float = cpl - cpv
ls1: float = xls + (cpi - cpv) * T00
ls2: float = cpi - cpv

kappa: float = (cp - cv) / cp
kappa_d: float = Rd / cp
rp00: float = 1.0 / p00
reps: float = Rv / Rd
eps: float = epsilon
rddcp: float = kappa_d
cpdrd: float = cp / Rd
cpdg: float = cp / g

converge: float = 0.0002
