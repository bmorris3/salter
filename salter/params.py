from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import batman
import numpy as np

from .cache import planet_props

__all__ = ['kic_to_params', 'transit_model']


def handle_w(w):
    if type(w) is np.ma.core.MaskedConstant:
        return 90.0
    else:
        return w

def handle_a(a):
    if type(a) is np.ma.core.MaskedConstant:
        return 90.0
    else:
        return a


def kic_to_params(kic):
    """
    Parameters
    ----------
    table :
    kic :

    Examples
    --------
    >>> from salter import kic_to_params
    >>> params = kic_to_params(9705459)
    """
    table = planet_props.table
    params = batman.TransitParams()       # object to store transit parameters
    params.t0 = table.loc[kic]['koi_time0bk'] + 2454833.0  # time of inferior conjunction
    params.per = table.loc[kic]['koi_period']                       # orbital period
    params.rp = np.sqrt(table.loc[kic]['koi_depth'] / 1e6)          # planet radius (in units of stellar radii)
    params.ecc = 0 #table.loc[kic]['koi_eccen']
    params.w = 90 #table.loc[kic]['koi_longp']

    params.duration = table.loc[kic]['koi_duration'] / 24.0  # [days]
    params.b = table.loc[kic]['koi_impact']  # impact parameter

    a_rs, inc = T14b2aRsi(params.per,params.duration, params.b,
                          params.rp, params.ecc, params.w)

    params.a = a_rs
    params.inc = inc

    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.2, 0.1]      #limb darkening coefficients
    return params


def transit_model(kic, times):
    """
    Parameters
    ----------
    table :
    kic :
    times :
    """
    table = planet_props.table
    params = kic_to_params(kic)
    m = batman.TransitModel(params, times)    #initializes model
    flux = m.light_curve(params)                    #calculates light curve
    return flux


def impact_parameter(transit_params):
    """
    Calculate impact parameter of transit from other transit parameters. From
    Winn 2010, Eqn 7 [1]_
    Parameters
    ----------
    transit_params : `~batman.TransitParams`
        Transit light curve parameters
    Returns
    -------
    b : float
        Impact parameter
    .. [1] http://adsabs.harvard.edu/abs/2010arXiv1001.2010W
    """
    e = transit_params.ecc  # eccentricity
    w = transit_params.w    # long. of pericenter [deg]
    a_on_Rs = transit_params.a  # a/R_s
    i = transit_params.inc  # inclination [deg]

    b = (a_on_Rs * np.cos(np.radians(i)) *
         (1 - e**2) / (1 + e*np.sin(np.radians(w))))
    return b


def T14b2aRsi(P, T14, b, RpRs, eccentricity, omega):
    '''
    Convert from duration and impact param to a/Rs and inclination
    '''
    beta = (1 - eccentricity**2)/(1 + eccentricity*np.sin(np.radians(omega)))
    C = np.sqrt(1 - eccentricity**2)/(1 + eccentricity*np.sin(np.radians(omega)))
    i = np.arctan(beta * np.sqrt((1 + RpRs)**2 - b**2)/(b*np.sin(T14*np.pi/(P*C))))
    aRs = b/(np.cos(i) * beta)
    return aRs, np.degrees(i)