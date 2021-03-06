from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os

import numpy as np
import kplr
from astropy.io import ascii
import h5py
from astropy.utils.console import ProgressBar
from astropy.utils.data import download_file
from astropy.table import Column, unique, join

__all__ = ['cache_light_curves', 'get_planets_table', 'cache_planets_table',
           'planet_props', 'lc_archive']


kic_numbers_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.path.pardir, 'data', 'kics.csv')

planet_table_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.path.pardir, 'data', 'joined_table.csv')

light_curves_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.path.pardir, 'data', 'light_curves.hdf5')

stats_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          os.path.pardir, 'data', 'stats.hdf5')


def cache_light_curves():
    """
    Run this after running `choose_targets.ipynb` in order to cache light curves
    into a local HDF5 archive.

    Examples
    --------
    >>> from salter import cache_light_curves; cache_light_curves()
    """
    if os.path.exists(light_curves_path):
        raise ValueError('Light curves file already exists, at {0}'
                         .format(light_curves_path))

    if not os.path.exists(kic_numbers_path):
        raise ValueError("You must first run the `choose_targets.ipynb` "
                         "notebook before running `cache_light_curves`")

    kics = ascii.read(kic_numbers_path, format='no_header')['col1']

    client = kplr.API()

    # Create archive
    f = h5py.File(light_curves_path, 'w')

    with ProgressBar(len(kics)) as bar:
        for kic in kics:
            if str(kic) not in f.keys():
                # Find a KIC
                star = client.star(kic)

                # Download the lightcurves for this KOI.
                lightcurves = star.get_light_curves(short_cadence=False)

                # Loop over the datasets and read in the data.
                time, flux, ferr, quality, quarter = [], [], [], [], []

                for i, lc in enumerate(lightcurves):
                    with lc.open() as lc_file:
                        # The lightcurve data are in the first FITS HDU.
                        hdu_data = lc_file[1].data

                        time.append(hdu_data["time"])
                        flux.append(hdu_data["sap_flux"])
                        ferr.append(hdu_data["sap_flux_err"])
                        quality.append(hdu_data["sap_quality"])
                        quarter.append(i * np.ones_like(hdu_data["time"]))

                data = np.vstack(list(map(np.concatenate, [time, flux, ferr, quality, quarter]))).T
                f.create_dataset(str(kic), data=data)
                f.flush()
                bar.update()
    f.close()


def cache_planets_table():
    """
    Cache a joined table containing data from the NASA Exoplanet Archive and
    the Exoplanet Orbit Database.

    To get the table, run the `~salter.get_planets_table()` function.
    """
    NEA_URL = 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative'
    EOD_URL = 'http://exoplanets.org/csv-files/exoplanets.csv'
    nea_table = ascii.read(download_file(NEA_URL, cache=False))
    eod_table = ascii.read(download_file(EOD_URL, cache=False))
    eod_table2 = eod_table[~eod_table['KEPID'].mask]
    nea_table2 = nea_table[~nea_table['kepid'].mask]

    eod_table2.add_column(Column(eod_table2['KEPID'], 'kepid'))

    joined_table = join(eod_table2, nea_table2, keys=['kepid'])

    ascii.write(joined_table, planet_table_path, format='csv')


def get_planets_table():
    """
    Get the joined planets table from the NASA Exoplanet Archive and
    the Exoplanet Orbit Database.

    Returns
    -------
    table : `~astropy.table.Table`
        Table of exoplanet properties
    """
    if not os.path.exists(planet_table_path):
        raise ValueError("You must run salter.cache.cache_planets_table first "
                         "before you can run get_joined_table")
    table = ascii.read(planet_table_path, format='csv')

    # Toss out multis
    first_kois_only = np.array([koi.endswith('01')
                                for koi in table['kepoi_name']])
    table = table[first_kois_only]
    table.add_index('kepid')

    # Ensure only unique results
    unique_table = unique(table, keys='kepid')
    unique_table.add_index('kepid')

    return unique_table


class PlanetProperties(object):
    """
    Cache manager for planet properties table.
    """
    def __init__(self):
        self._table = None

    @property
    def table(self):
        """
        Column definitions can be found at [1]_ and [2]_.

        References
        ----------
        .. [1] http://exoplanets.org/help/common/data
        .. [2] https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html
        """
        if self._table is None:
            self._table = get_planets_table()
        return self._table


class LightCurveArchive(object):
    """
    Light curve HDF5 archive manager
    """
    def __init__(self):
        self._file = None

    @property
    def file(self):
        """
        Return an open HDF5 file stream of the light curve archive.
        """
        if self._file is None:
            self._file = h5py.File(light_curves_path, 'r')
        return self._file


planet_props = PlanetProperties()
lc_archive = LightCurveArchive()
