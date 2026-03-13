#!/usr/bin/env python
"""This module provides dynamics calculations for atmospheric analysis."""

from __future__ import annotations

import math
import numpy as np
import numpy.typing as npt
from pymeteo.constants import gravity, km2m, missingval
import pymeteo.interp

ArrayLike = npt.NDArray[np.floating]


# Mean wind in the layer [zmin,zmax]
def avg_wind(
    _u: ArrayLike, _v: ArrayLike, _z: ArrayLike, zmin: float, zmax: float
) -> tuple[float, float]:
    """Calculate mean wind in a layer.

    :param _u: U wind component array
    :param _v: V wind component array
    :param _z: Height array (m)
    :param zmin: Bottom of layer (m)
    :param zmax: Top of layer (m)
    :returns: Tuple of (mean u, mean v); (0.0, 0.0) if no data in layer
    """
    # TODO: np.average !
    # REPLACED BY mean_wind below
    u = 0.0
    v = 0.0
    n = 0
    for i in np.arange(0, len(_z), 1):
        if _z[i] > zmin and _z[i] < zmax:
            u += _u[i]
            v += _v[i]
            n += 1
    if n == 0:
        return (0.0, 0.0)
    return (u / n, v / n)


def storm_motion_rasmussen(
    _u: ArrayLike, _v: ArrayLike, _z: ArrayLike
) -> tuple[float, float, float, float]:
    """Calculate storm motion using Rasmussen (1984) method.

    :param _u: U wind component array
    :param _v: V wind component array
    :param _z: Height array (m)
    :returns: Tuple of (u_right, v_right, u_left, v_left) storm motion components
    """
    # rasmussen (1984)
    u_0_500 = avg_wind(_u, _v, _z, 0.0, 500.0)
    u_4km = avg_wind(_u, _v, _z, 3500.0, 4500.0)

    dist_60pct = (
        math.sqrt((u_4km[0] - u_0_500[0]) ** 2 + (u_4km[1] - u_0_500[1]) ** 2) * 0.6
    )
    du = u_4km[0] - u_0_500[0]
    dv = u_4km[1] - u_0_500[1]
    theta = math.atan2(dv, du)
    u60 = dist_60pct * math.cos(theta) + u_0_500[0]
    v60 = dist_60pct * math.sin(theta) + u_0_500[1]
    theta -= math.pi / 2.0
    dist = 8.7
    u_cr = dist * math.cos(theta) + u60
    v_cr = dist * math.sin(theta) + v60
    theta += math.pi
    u_cl = dist * math.cos(theta) + u60
    v_cl = dist * math.sin(theta) + v60

    return u_cr, v_cr, u_cl, v_cl


def storm_motion_bunkers(
    _u: ArrayLike, _v: ArrayLike, _z: ArrayLike
) -> tuple[float, float, float, float]:
    """Calculate storm motion using Bunkers (2000) method.

    :param _u: U wind component array
    :param _v: V wind component array
    :param _z: Height array (m)
    :returns: Tuple of (u_right, v_right, u_left, v_left) storm motion components
    """
    # bunkers (2000)
    u_0_500 = avg_wind(_u, _v, _z, 0.0, 500.0)
    u_0_6km = avg_wind(_u, _v, _z, 0.0, 6000.0)
    du = u_0_6km[0] - u_0_500[0]
    dv = u_0_6km[1] - u_0_500[1]
    theta = math.atan2(dv, du) - math.pi / 2.0
    dist = 7.5
    u_cr = dist * math.cos(theta) + u_0_6km[0]
    v_cr = dist * math.sin(theta) + u_0_6km[1]
    theta += math.pi
    u_cl = dist * math.cos(theta) + u_0_6km[0]
    v_cl = dist * math.sin(theta) + u_0_6km[1]

    return u_cr, v_cr, u_cl, v_cl


def circulation(
    u: ArrayLike,
    v: ArrayLike,
    w: ArrayLike,
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
) -> tuple[float, ArrayLike, ArrayLike]:
    """Calculate circulation around a closed circuit.

    :param u: U wind component array
    :param v: V wind component array
    :param w: W wind component array
    :param x: X position array (km)
    :param y: Y position array (km)
    :param z: Z position array (km)
    :returns: Tuple of (circulation, vdotdl, normalized vdotdl)
    """
    n = len(x)

    # local variables
    uavg = np.zeros((n), np.float32)
    vavg = np.zeros((n), np.float32)
    wavg = np.zeros((n), np.float32)

    # initialize returned variables
    vdotdl1 = np.empty((n), np.float32)
    vdotdl2 = np.empty((n), np.float32)
    C = 0.0

    if (
        (u == missingval).any()
        or (v == missingval).any()
        or (w == missingval).any()
        or (x == missingval).any()
        or (y == missingval).any()
        or (z == missingval).any()
    ):
        del uavg, vavg, wavg, vdotdl1, vdotdl2
        return (missingval, np.zeros(0, np.float32), np.zeros(0, np.float32))

    # calculate wind and dl around the circuit
    uavg[0 : n - 1] = 0.5 * (u[0 : n - 1] + u[1:n])
    vavg[0 : n - 1] = 0.5 * (v[0 : n - 1] + v[1:n])
    wavg[0 : n - 1] = 0.5 * (w[0 : n - 1] + w[1:n])

    uavg[n - 1] = 0.5 * (u[0] + u[n - 1])
    vavg[n - 1] = 0.5 * (v[0] + v[n - 1])
    wavg[n - 1] = 0.5 * (w[0] + w[n - 1])

    # ediff1d returns an array of the difference between
    # each element of the passed array, which is to say,
    # it returns the delta of each array element.  perfect.
    dx = np.ediff1d(x, to_end=x[0] - x[n - 1])
    dy = np.ediff1d(y, to_end=y[0] - y[n - 1])
    dz = np.ediff1d(z, to_end=z[0] - z[n - 1])

    # assumes clockwise parcels
    C = np.sum(-uavg * dx - vavg * dy - wavg * dz)
    vdotdl1 = -uavg * dx - vavg * dy - wavg * dz
    vdotdl2 = vdotdl1 / (dx**2 + dy**2 + dz**2) ** 0.5
    vdotdl1 = vdotdl1 * km2m

    C = C * km2m  # C now in units m2/s

    del uavg, vavg, wavg, dx, dy, dz
    return (C, vdotdl1, vdotdl2)


def integral_Bdz(th: ArrayLike, thp: ArrayLike, z: ArrayLike) -> float:
    """Calculate the integral of buoyancy over height.

    :param th: Environmental potential temperature array (K)
    :param thp: Parcel potential temperature array (K)
    :param z: Height array (km)
    :returns: Integral of buoyancy * dz (m2/s2)
    """
    n = len(z)

    th_avg = np.empty((n), np.float32)
    thp_avg = np.empty((n), np.float32)
    dz = np.empty((n), np.float32)

    intBdz = 0.0

    if (th == missingval).any() or (thp == missingval).any() or (z == missingval).any():
        intBdz = missingval
        del th_avg, thp_avg, dz
        return intBdz

    th_avg[0 : n - 1] = 0.5 * (th[0 : n - 1] + th[1:n])
    thp_avg[0 : n - 1] = 0.5 * (thp[0 : n - 1] + thp[1:n])

    th_avg[n - 1] = 0.5 * (th[0] + th[n - 1])
    thp_avg[n - 1] = 0.5 * (thp[0] + thp[n - 1])

    dz = np.ediff1d(z, to_end=z[0] - z[n - 1])

    intBdz = np.sum(-dz * gravity * thp_avg / (th_avg - thp_avg))

    del th_avg, thp_avg, dz
    intBdz = intBdz * km2m  # units now m2/s2

    return intBdz


def integral_dt(i: ArrayLike, t: ArrayLike) -> float:
    """Calculate the integral of a quantity over time.

    :param i: Integrand array
    :param t: Time array
    :returns: Integral value
    """
    n = len(t)

    iavg = np.empty((n - 1), np.float32)
    dt = np.empty((n - 1), np.float32)

    iavg[0 : n - 1] = 0.5 * (i[0 : n - 1] + i[1:n])

    dt = np.ediff1d(t)

    integral = np.sum(iavg * dt)

    return integral


# helper functions


def uv_to_deg(u: float, v: float) -> tuple[float, float]:
    """Transform u, v wind components to direction and magnitude.

    :param u: U wind component (m/s)
    :param v: V wind component (m/s)
    :returns: Tuple of (wind direction in degrees, wind speed in m/s)
    """
    direction = np.arctan2(u, v) * (180.0 / np.pi)
    speed = (u**2 + v**2) ** (0.5)

    return direction, speed


def wind_deg_to_uv(direction: float, speed: float) -> tuple[float, float]:
    """Convert direction and speed into u, v wind components.

    :param direction: Wind direction (mathematical angle, degrees)
    :param speed: Wind magnitude (m/s)
    :returns: Tuple of (u, v) wind components (m/s)
    """
    u = speed * np.sin(np.pi * (direction + 180.0) / 180.0)
    v = speed * np.cos(np.pi * (direction + 180.0) / 180.0)

    return u, v


def shear(
    _u: ArrayLike, _v: ArrayLike, _z: ArrayLike, zbot: float, ztop: float
) -> tuple[float, float]:
    """Calculate the wind shear in a layer.

    :param _u: U wind component array (m/s)
    :param _v: V wind component array (m/s)
    :param _z: Height array (m)
    :param zbot: Bottom of the layer (m)
    :param ztop: Top of the layer (m)
    :returns: Tuple of (u_shear, v_shear)
    """
    if zbot < _z[0]:
        zbot = _z[0]

    ubot = pymeteo.interp.linear(_z, _u, zbot)
    vbot = pymeteo.interp.linear(_z, _v, zbot)
    utop = pymeteo.interp.linear(_z, _u, ztop)
    vtop = pymeteo.interp.linear(_z, _v, ztop)

    u = utop - ubot
    v = vtop - vbot
    return u, v


def srh(
    _u: ArrayLike,
    _v: ArrayLike,
    _z: ArrayLike,
    zbot: float,
    ztop: float,
    cx: float,
    cy: float,
) -> float:
    """Calculate storm relative helicity in a layer.

    :param _u: U wind component array (m/s)
    :param _v: V wind component array (m/s)
    :param _z: Height array (m)
    :param zbot: Bottom of the layer (m)
    :param ztop: Top of the layer (m)
    :param cx: U component of storm motion (m/s)
    :param cy: V component of storm motion (m/s)
    :returns: Storm relative helicity (m2/s2)
    """
    if zbot < _z[0]:
        zbot = _z[0]

    dz = 10.0
    z = np.arange(zbot, ztop + dz, dz)
    nk = len(z)
    u = np.empty(nk, np.float32)
    v = np.empty(nk, np.float32)
    du = np.empty(nk - 1, np.float32)
    dv = np.empty(nk - 1, np.float32)
    uavg = np.empty(nk - 1, np.float32)
    vavg = np.empty(nk - 1, np.float32)

    for k in range(nk):
        u[k] = pymeteo.interp.linear(_z, _u, z[k])
        v[k] = pymeteo.interp.linear(_z, _v, z[k])

    du = np.ediff1d(u)
    dv = np.ediff1d(v)

    uavg[0 : nk - 1] = 0.5 * (u[0 : nk - 1] + u[1:nk])
    vavg[0 : nk - 1] = 0.5 * (v[0 : nk - 1] + v[1:nk])

    srh = np.sum(-(uavg - cx) * dv + (vavg - cy) * du)
    return srh


def mean_wind(
    _u: ArrayLike, _v: ArrayLike, _z: ArrayLike, zbot: float, ztop: float
) -> tuple[float, float]:
    """Calculate the mean wind in a layer.

    :param _u: U wind component array (m/s)
    :param _v: V wind component array (m/s)
    :param _z: Height array (m)
    :param zbot: Bottom of the layer (m)
    :param ztop: Top of the layer (m)
    :returns: Tuple of (mean u, mean v) (m/s)
    """
    if zbot < _z[0]:
        zbot = _z[0]

    dz = 10.0
    z = np.arange(zbot, ztop + dz, dz)
    nk = len(z)
    u = np.empty(nk, np.float32)
    v = np.empty(nk, np.float32)

    for k in range(nk):
        u[k] = pymeteo.interp.linear(_z, _u, z[k])
        v[k] = pymeteo.interp.linear(_z, _v, z[k])

    uavg = np.mean(u, dtype=np.float64)
    vavg = np.mean(v, dtype=np.float64)
    return uavg, vavg


def brn(_u: ArrayLike, _v: ArrayLike, _z: ArrayLike, cape: float) -> float:
    """Calculate the Bulk Richardson Number.

    :param _u: U wind component array (m/s)
    :param _v: V wind component array (m/s)
    :param _z: Height array (m)
    :param cape: Convective Available Potential Energy (J/kg)
    :returns: Bulk Richardson Number (dimensionless)
    """
    u06avg = mean_wind(_u, _v, _z, 0.0, 6000.0)
    u0500avg = mean_wind(_u, _v, _z, 0.0, 500.0)
    u = u06avg[0] - u0500avg[0]
    v = u06avg[1] - u0500avg[1]

    brn = cape / (0.5 * (u**2 + v**2))
    return brn
