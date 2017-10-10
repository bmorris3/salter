from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import batman

__all__ = ['kic_to_params']


def kic_to_params(table, kic):
    """
    Parameters
    ----------
    table :
    kic :

    Examples
    --------
    >>> from salter import get_joined_table, kic_to_params
    >>> table = get_joined_table()
    >>> params = kic_to_params(table, 9705459)
    """
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = table.loc[kic]['koi_time0bk'] + 2454833.0                        #time of inferior conjunction
    params.per = table.loc[kic]['koi_period']                       #orbital period
    params.rp = table.loc[kic]['RR']                      #planet radius (in units of stellar radii)
    params.a = table.loc[kic]['A']                        #semi-major axis (in units of stellar radii)
    params.inc = table.loc[kic]['I']                      #orbital inclination (in degrees)
    params.ecc = table.loc[kic]['ECC']                       #eccentricity
    params.w = table.loc[kic]['OM']                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.2, 0.1]      #limb darkening coefficients
    return params


