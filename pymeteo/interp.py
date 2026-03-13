"""Interpolation utilities for atmospheric profile data."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pymeteo.constants

# Type alias for array-like inputs
ArrayLike = npt.NDArray[np.floating]


def linear(dim: ArrayLike, var: ArrayLike, dval: float) -> float:
    """Interpolate var to the value dval along dimension dim.

    :param dim: Dimension array (e.g. height)
    :param var: Variable array to interpolate
    :param dval: Target dimension value to interpolate to
    :returns: Interpolated value of var at dval
    """
    if len(dim) != len(var):
        raise Exception("Dimensions of dim and var do not match!")
    dim = np.ma.masked_invalid(dim)
    var = np.ma.masked_invalid(var)
    array_mask = dim.mask | var.mask
    array_x = dim.data[~array_mask]
    array_y = var.data[~array_mask]
    index_sort = np.argsort(array_x)

    return np.interp(dval, array_x[index_sort], array_y[index_sort])


def interp_height(z: ArrayLike, p: ArrayLike, plvl: float) -> float | int:
    """Interpolate height to a pressure level.

    :param z: Height array (m)
    :param p: Pressure array (Pa)
    :param plvl: Target pressure level (Pa)
    :returns: Interpolated height (m), 0 if below surface, -1 if above top
    """
    nlevs = len(p)

    # check bounds
    if plvl > p[0]:
        return 0
    if plvl < p[nlevs - 1]:
        return -1

    z0 = 0
    while p[z0] > plvl:
        z0 = z0 + 1

    if p[z0] == plvl:
        return z[z0]

    z1 = nlevs - 1
    while p[z1] < plvl:
        z1 = z1 - 1

    # interpolate to height.
    # Code adapted from NSHARP 95 John Hart NSSFC KCMO.

    zdiff = z[z1] - z[z0]
    pdiff = np.log(p[z0] / p[z1])
    pdist = np.log(p[z0] / plvl)
    height = z[z0] + ((pdist / pdiff) * zdiff)

    return height


def interp_pressure(p: ArrayLike, z: ArrayLike, zlvl: float) -> float | None:
    """Interpolate pressure to a height level.

    :param p: Pressure array (Pa)
    :param z: Height array (m)
    :param zlvl: Target height level (m)
    :returns: Interpolated pressure (Pa), or None if above top
    """
    nlevs = len(p)

    z0 = nlevs - 1
    while z[z0] > zlvl:
        z0 = z0 - 1

    if z[z0] == zlvl:
        return p[z0]

    z1 = 0
    while z[z1] < zlvl:
        z1 = z1 + 1
        if z1 >= nlevs:
            return None

    zdiff = z[z1] - z[z0]
    zdist = zlvl - z[z0]
    pdiff = np.log(p[z1] / p[z0])
    pres = p[z0] * np.exp((zdist / zdiff) * pdiff)
    return pres
