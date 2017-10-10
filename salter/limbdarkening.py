# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astroquery.vizier import Vizier

__all__ = ['quad', 'ld_table']


class LDTable(object):
    def __init__(self):
        self._table = None

    @property
    def table(self):
        if self._table is None:

            table = Vizier.get_catalogs('J/A+A/552/A16/limb1-4')[0]
            self._table = table
        return self._table

ld_table = LDTable()


def get_closest_model(Teff, logg, filt):
    '''
    All grid points have metallicities Z=0 (solar), microturbulence xi=2

    Parameters
    ----------
    filt : string
        Name of the filter used.

    logg : float
        Logarithm of the stellar surface gravity.

    Teff : float
        Host star effective temperature.
    '''
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
        closestmodel = get_closest_model(round(Teff, -2), logg, filt)
    else:
        closestmodel = get_closest_model(5800, logg, filt)

    return ld_table.table['a'][closestmodel], ld_table.table['b'][closestmodel]
