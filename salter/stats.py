from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from scipy.stats import ttest_ind, anderson_ksamp, ks_2samp
from .lightcurve import concatenate_transit_light_curves
import numpy as np
import matplotlib.pyplot as plt

__all__ = ['Residuals']


class Residuals(object):
    def __init__(self, transits, params, buffer_duration=0.15):
        """
        Transit light curve residuals.

        Parameters
        ----------
        transits : list of (or single) `~salter.TransitLightCurve` objects
            list of transits
        params : `~batman.TransitParams()`
            transiting planet parameters
        buffer_duration : float
            fraction of transit duration to ignore centered on ingress and
            egress.
        """
        if type(transits) is list:
            all_transits = concatenate_transit_light_curves(transits)
        else:
            all_transits = transits

        self.residuals = all_transits.fluxes - all_transits.transit_model()
        self.params = params
        self.egress_phase = params.duration/2/params.per
        self.phases = all_transits.phases()

        mask_nans = np.isnan(self.residuals)

        buffer = buffer_duration * params.duration / params.per

        out_of_transit = (((-self.egress_phase - buffer/2 > all_transits.phases()) |
                           (self.egress_phase + buffer/2 < all_transits.phases())) &
                          np.logical_not(mask_nans))
        in_transit = (((-self.egress_phase + buffer/2 < all_transits.phases()) &
                           (self.egress_phase - buffer/2 > all_transits.phases())) &
                          np.logical_not(mask_nans))

        self.out_of_transit = out_of_transit
        self.in_transit = in_transit
        self.before_midtransit = self.phases < 0
        self.after_midtransit = self.phases > 0

    def plot(self):
        """
        Generate a quick plot of the transit residuals.

        Returns
        -------
        fig : `~matplotlib.pyplot.Figure`
            Figure object
        ax : `~matplotlib.pyplot.Axes`
            axis object
        """
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
        Independent two-sample T test from `~scipy.stats.ttest_ind`.

        Parameters
        ----------
        attrs1 : list of attributes
            List of conditions in first sample
        attrs2 : list of attributes
            List of conditions in second sample

        Returns
        -------
        pvalue : float
            p value.

        Examples
        --------
        >>> import numpy as np
        >>> import batman
        >>> from salter import LightCurve
        >>> # Create example transiting planet properties
        >>> params = batman.TransitParams()
        >>> params.t0 = 0.5
        >>> params.rp = 0.1
        >>> params.per = 1
        >>> params.duration = 0.3
        >>> params.inc = 90
        >>> params.w = 90
        >>> params.ecc = 0
        >>> params.a = 10
        >>> params.limb_dark = 'quadratic'
        >>> params.u = [0.2, 0.1]
        >>> # Create example transit light curves:
        >>> transits = [LightCurve(times=i + np.linspace(0, 1, 500),
        >>>                        fluxes=np.random.randn(500),
        >>>                        params=params) for i in range(10)]
        >>> # Create residuals object
        >>> r = Residuals(transits, params)
        >>> # How significant is the difference between the means of the fluxes in and out-of-transit?
        >>> r.ttest('out_of_transit', 'in_transit')
        0.310504218041
        >>> # How significant is the difference between the means of the in-transit fluxes before and after midtransit?
        >>> r.ttest(['in_transit', 'before_midtransit'], ['in_transit', 'after_midtransit'])
        0.823997471194
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        return ttest_ind(sample1, sample2, equal_var=False).pvalue

    def ks(self, attrs1, attrs2):
        """
        Two-sample KS test from `~scipy.stats.ks_2samp`

        Parameters
        ----------
        attrs1 : list of attributes
            List of conditions in first sample
        attrs2 : list of attributes
            List of conditions in second sample

        Returns
        -------
        pvalue : float
            p value.

        Examples
        --------
        >>> import numpy as np
        >>> import batman
        >>> from salter import LightCurve
        >>> # Create example transiting planet properties
        >>> params = batman.TransitParams()
        >>> params.t0 = 0.5
        >>> params.rp = 0.1
        >>> params.per = 1
        >>> params.duration = 0.3
        >>> params.inc = 90
        >>> params.w = 90
        >>> params.ecc = 0
        >>> params.a = 10
        >>> params.limb_dark = 'quadratic'
        >>> params.u = [0.2, 0.1]
        >>> # Create example transit light curves:
        >>> transits = [LightCurve(times=i + np.linspace(0, 1, 500),
        >>>                        fluxes=np.random.randn(500),
        >>>                        params=params) for i in range(10)]
        >>> r = Residuals(transits, params)
        >>> # How significant is the difference between the distributions of the fluxes in and out-of-transit?
        >>> r.ks('out_of_transit', 'in_transit')
        0.91710727901331124
        >>> # How significant is the difference between the distributions of the in-transit fluxes before and after midtransit?
        >>> r.ks(['in_transit', 'before_midtransit'], ['in_transit', 'after_midtransit'])
        0.39171715554793468
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        return ks_2samp(sample1, sample2).pvalue

    def anderson(self, attrs1, attrs2):
        """
        k-sample Anderson test from `~scipy.stats.anderson_ksamp`.

        Parameters
        ----------
        attrs1 : list of attributes
            List of conditions in first sample
        attrs2 : list of attributes
            List of conditions in second sample

        Returns
        -------
        sig : float
            significance level (see `~scipy.stats.anderson_ksamp`)

        Examples
        --------
        >>> import numpy as np
        >>> import batman
        >>> from salter import LightCurve
        >>> # Create example transiting planet properties
        >>> params = batman.TransitParams()
        >>> params.t0 = 0.5
        >>> params.rp = 0.1
        >>> params.per = 1
        >>> params.duration = 0.3
        >>> params.inc = 90
        >>> params.w = 90
        >>> params.ecc = 0
        >>> params.a = 10
        >>> params.limb_dark = 'quadratic'
        >>> params.u = [0.2, 0.1]
        >>> # Create example transit light curves:
        >>> transits = [LightCurve(times=i + np.linspace(0, 1, 500),
        >>>                        fluxes=np.random.randn(500),
        >>>                        params=params) for i in range(10)]
        >>> r = Residuals(transits, params)
        >>> # How significant is the difference between the distributions of the fluxes in and out-of-transit?
        >>> r.anderson('out_of_transit', 'in_transit')
        1.1428634099527666
        >>> # How significant is the difference between the distributions of the in-transit fluxes before and after midtransit?
        >>> r.anderson(['in_transit', 'before_midtransit'], ['in_transit', 'after_midtransit'])
        0.2792395871784852
        """
        sample1, sample2 = self._and_reduce(attrs1, attrs2)

        try:
            return anderson_ksamp([sample1, sample2]).significance_level
        except OverflowError:
            return 0
