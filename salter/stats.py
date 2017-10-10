from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from scipy.stats import ttest_ind, anderson_ksamp, ks_2samp
from .lightcurve import concatenate_transit_light_curves
import numpy as np
import matplotlib.pyplot as plt

__all__ = ['Residuals']


class Residuals(object):
    def __init__(self, transits, params):
        """
        Residuals container object

        Parameters
        ----------
        transits : list
        params : `~batman.TransitParams()`
        """

        all_transits = concatenate_transit_light_curves(transits)

        self.residuals = all_transits.fluxes - all_transits.transit_model()
        self.params = params
        self.egress_phase = params.duration/2/params.per
        self.phases = all_transits.phases()

        mask_nans = np.isnan(self.residuals)

        out_of_transit = (((-self.egress_phase > all_transits.phases()) |
                           (self.egress_phase < all_transits.phases())) &
                          np.logical_not(mask_nans))
        in_transit = np.logical_not(out_of_transit) & np.logical_not(mask_nans)

        self.out_of_transit = out_of_transit
        self.in_transit = in_transit
        self.before_midtransit = self.phases < 0
        self.after_midtransit = self.phases > 0

    def plot(self):
        fig, ax = plt.subplots()

        ax.plot(self.phases[self.out_of_transit], self.residuals[self.out_of_transit], 'k.', alpha=0.3)
        ax.plot(self.phases[self.in_transit], self.residuals[self.in_transit], 'r.', alpha=0.3)
        return fig, ax

    def _and_reduce(self, attrs1, attrs2):
        """
        If ``attrs`` is a list of attributes, AND-combine all attributes.
        If ``attrs`` is a single attribute, then just return that one.
        """
        if type(attrs1) is list:
            condition1 = np.logical_and.reduce([getattr(self, attr)
                                             for attr in attrs1])
        else:
            condition1 = getattr(self, attrs1)

        if type(attrs2) is list:
            condition2 = np.logical_and.reduce([getattr(self, attr)
                                             for attr in attrs2])
        else:
            condition2 = getattr(self, attrs2)

        return self.residuals[condition1], self.residuals[condition2]

    def ttest(self, attrs1, attrs2):
        """
        Independent two-sample T test.

        Parameters
        ----------
        attrs1 : list of attributes
        attrs2 : list of attributes
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        return ttest_ind(sample1, sample2, equal_var=False).pvalue

    def ks(self, attrs1, attrs2):
        """
        Two-sample KS test.

        Parameters
        ----------
        attrs1 : list of attributes
        attrs2 : list of attributes
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        return ks_2samp(sample1, sample2).pvalue

    def anderson(self, attrs1, attrs2):
        """
        k-sample Anderson test

        Parameters
        ----------
        attrs1 : list of attributes
        attrs2 : list of attributes
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        return anderson_ksamp([sample1, sample2]).significance_level
