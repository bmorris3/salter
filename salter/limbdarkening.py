# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astroquery.vizier import Vizier

__all__ = ['quad', 'ld_table', 'q2u', 'u2q']


class LDTable(object):
    """
    Limb darkening table object.
    """
    def __init__(self):
        self._table = None

    @property
    def table(self):
        if self._table is None:

            table = Vizier.get_catalogs('J/A+A/552/A16/limb1-4')[0]
            self._table = table
        return self._table

ld_table = LDTable()


def _get_closest_model(Teff, logg, filt):
    """
    All grid points have metallicities Z=0 (solar), microturbulence xi=2

    Parameters
    ----------
    filt : string
        Name of the filter used.

    logg : float
        Logarithm of the stellar surface gravity.

    Teff : float
        Host star effective temperature.
    """
    correct_filter = ld_table.table['Filt'] == filt
    method = ld_table.table['Met'] == 'F'
    mask = (correct_filter & method)
    return np.argmin(np.abs(ld_table.table['Teff'] - Teff)/1000 +
                     np.abs(ld_table.table['logg'] - logg) +
                     (~mask)*1000)


def quad(Teff, logg, filt):
    """
    Quadratic limb-darkening law.
    All grid points have metallicities Z=0 (solar), microturbulence xi=2.

    Parameters
    ----------
    filt : string
        Name of the filter used.

    logg : float
        Logarithm of the stellar surface gravity.

    Teff : float
        Host star effective temperature.

    method : string, optional. Either 'F' or 'L'
        Method of computation: least-square or flux conservation

    Returns
    -------
    a, b : float, float
        The linear and quadratic limb-darkening terms for a quadratic law:
        I(mu)/I(1) = 1 - a*(1 - mu) - b*(1 - mu)**2
    """
    if type(Teff) is int or type(Teff) is float:
        closestmodel = _get_closest_model(round(Teff, -2), logg, filt)
    else:
        closestmodel = _get_closest_model(5800, logg, filt)

    return ld_table.table['a'][closestmodel], ld_table.table['b'][closestmodel]


def u2q(u1, u2, warnings=True):
    """
    Convert the linear and quadratic terms of the quadratic limb-darkening
    parameterization -- called `u_1` and `u_2` in Kipping 2013 or `a` and `b` in
    Claret et al. 2013 -- and convert them to `q_1` and `q_2` as described in
    Kipping 2013:

    http://adsabs.harvard.edu/abs/2013MNRAS.435.2152K

    Parameters
    ----------
    u1 : float
        Linear component of quadratic limb-darkening

    u2 : float
        Quadratic component of quadratic limb-darkening

    Returns
    -------
    (q1, q2) : tuple of floats
        Kipping (2013) style quadratic limb-darkening parameters
    """
    q1 = (u1 + u2)**2
    q2 = 0.5*u1/(u1+u2)
    if warnings and (u1 < 0 or u2 < 0):
        print("WARNING: The quadratic limb-darkening parameters " +
              "u1={0:.3f} or u2={0:.3f} violate Kipping's ".format(u1, u2) +
              "conditions for a monotonically increasing or everywhere-" +
              "positive intensity profile. Returning them as is.")
    return q1, q2


def q2u(q1, q2):
    """
    Convert the two parameter quadratic terms of the Kipping 2013 limb-
    darkening parameterization `q_1` and `q_2` to the standard linear and
    quadratic terms of the quadratic limb-darkening parameterization of
    Claret et al. 2013 -- called `u_1` and `u_2` in Kipping 2013 or `a` and `b` in
    Claret et al. 2013:

    http://adsabs.harvard.edu/abs/2013A%26A...552A..16C

    Parameters
    ----------
    q1 : float
        First component of Kipping 2013 quadratic limb-darkening

    q2 : float
        Second component of Kipping 2013 quadratic limb-darkening

    Returns
    -------
    (u1, u2) : tuple of floats
        Claret et al. (2013) style quadratic limb-darkening parameters
    """
    u1 = 2*np.sqrt(q1)*q2
    u2 = np.sqrt(q1)*(1-2*q2)
    return u1, u2
