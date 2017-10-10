from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import numpy as np
import kplr
from astropy.io import ascii
import os
import h5py
from astropy.utils.console import ProgressBar
from astropy.utils.data import download_file
from astropy.table import join
from astropy.table import Column

__all__ = ['cache_light_curves', 'get_planets_table', 'cache_joined_table']

kic_numbers_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.path.pardir, 'data', 'kics.csv')

planet_table_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.path.pardir, 'data', 'joined_table.csv')

light_curves_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.path.pardir, 'data', 'light_curves.hdf5')


def cache_light_curves():
    """
    Run:
    python -c "from salter import cache_light_curves; cache_light_curves()"
    """
    if os.path.exists(light_curves_path):
        raise ValueError('Light curves file already exists, at {0}'
                         .format(light_curves_path))

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
                time, flux, ferr, quality = [], [], [], []
                n_columns = 4

                for lc in lightcurves:
                    with lc.open() as lc_file:
                        # The lightcurve data are in the first FITS HDU.
                        hdu_data = lc_file[1].data

                        time.append(hdu_data["time"])
                        flux.append(hdu_data["sap_flux"])
                        ferr.append(hdu_data["sap_flux_err"])
                        quality.append(hdu_data["sap_quality"])

                data = np.vstack(list(map(np.concatenate, [time, flux, ferr, quality]))).T
                f.create_dataset(str(kic), data=data)
                f.flush()
                bar.update()


def cache_joined_table():
    """
    Cache a joined table containing data from the NASA Exoplanet Archive and
    the Exoplanet Orbit Database.

    To get the table, run the ``get_joined_table()`` function
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
    if not os.path.exists(planet_table_path):
        raise ValueError("You must run salter.cache.cache_joined_table first "
                         "before you can run get_joined_table")
    table = ascii.read(planet_table_path, format='csv')

    # Toss out multis
    first_kois_only = np.array([koi.endswith('01')
                                for koi in table['kepoi_name']])
    table = table[first_kois_only]
    table.add_index('kepid')
    return table
