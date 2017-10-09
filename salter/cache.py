from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import numpy as np
import kplr
from astropy.io import ascii
import os
import h5py
from astropy.utils.console import ProgressBar

__all__ = ['get_light_curves']

kic_numbers_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.path.pardir, 'data', 'kics.csv')


def get_light_curves():
    kics = ascii.read(kic_numbers_path, format='no_header')['col1']

    client = kplr.API()

    # Create archive
    f = h5py.File('data/light_curves.hdf5', 'w')

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