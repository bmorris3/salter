# Licensed under the MIT License - see LICENSE.rst
"""
Methods for taking the raw light curves from MAST and producing cleaned light
curves.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from copy import deepcopy

from astropy.time import Time
import astropy.units as u
import os
import numpy as np
import matplotlib.pyplot as plt
import batman
from scipy.ndimage import gaussian_filter
from scipy.optimize import fmin_l_bfgs_b

from .params import kic_to_params
from .limbdarkening import q2u, u2q
from .cache import lc_archive

__all__ = ['LightCurve', 'concatenate_transit_light_curves',
           'TransitLightCurve', 'concatenate_light_curves',
           'subtract_add_divide']

LCARCHIVE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              os.path.pardir, 'data', 'light_curves.hdf5')


def generate_lc_depth(times, depth, transit_params, exp_time=30/60/24):
    """
    Generate a model transit light curve.

    Parameters
    ----------
    times : `~numpy.ndarray`
        Times in JD
    depth : float
        Set depth independently from the setting in `transit_params`
    transit_params : `~batman.TransitParams`
        Transit light curve parameters

    Returns
    -------

    """
    #exp_time = (30*u.min).to(u.day).value

    transit_params.rp = np.sqrt(depth)

    m = batman.TransitModel(transit_params, times, supersample_factor=7,
                            exp_time=exp_time)
    model_flux = m.light_curve(transit_params)
    return model_flux


def lc_3param(times, rp, q1, q2, transit_params, exp_time=30/60/24):
    """
    Generate a model transit light curve.

    Parameters
    ----------
    times : `~numpy.ndarray`
        Times in JD
    rp : float
        Planet-to-star radius ratio
    q1 : float on interval [0, 1]
        Limb darkening parameter from Kipping 2013
    q2 : float on interval [0, 1]
        Limb darkening parameter from Kipping 2013
    transit_params : `~batman.TransitParams`
        Transit light curve parameters

    Returns
    -------
    flux : `~numpy.ndarray`
        Model light curve at ``times``
    """
    params = deepcopy(transit_params)

    params.rp = rp
    u1, u2 = q2u(q1, q2)
    params.u = [u1, u2]

    m = batman.TransitModel(params, times, supersample_factor=7,
                            exp_time=exp_time)
    model_flux = m.light_curve(params)
    return model_flux



class LightCurve(object):
    """
    Container object for light curves.
    """
    def __init__(self, times=None, fluxes=None, errors=None, quarters=None,
                 name=None, params=None):
        """
        Parameters
        ----------
        times : `~numpy.ndarray`
            Times in JD
        fluxes : `~numpy.ndarray`
            Fluxes (normalized or not)
        errors : `~numpy.ndarray`
            Uncertainties on the fluxes
        quarters : `~numpy.ndarray` (optional)
            Kepler Quarter for each flux
        name : str
            Name this light curve (optional)
        params : `~batman.TransitParams`
            Planet transit parameters
        """
        # if len(times) < 1:
        #    raise ValueError("Input `times` have no length.")

        if isinstance(times[0], Time) and isinstance(times, np.ndarray):
            times = Time(times)
        elif not isinstance(times, Time):
            times = Time(times, format='jd')

        self.times = times
        self.fluxes = fluxes
        if self.times is not None and errors is None:
            errors = np.zeros_like(self.fluxes) - 1
        self.errors = errors
        if self.times is not None and quarters is None:
            quarters = np.zeros_like(self.fluxes) - 1
        self.quarters = quarters
        self.name = name if type(name) is str else str(name)
        self.params = params

    @classmethod
    def from_hdf5(cls, kic):
        """
        Load a light curve from the HDF5 light curve archive with KIC number
        ``kic``.

        Parameters
        ----------
        kic: float
            KIC number
        """
        hdf5_file = lc_archive.file
        mask_nans = np.logical_not(np.isnan(hdf5_file[str(kic)][:, 0]))

        name = str(kic)

        params = kic_to_params(kic)

        return cls(times=hdf5_file[str(kic)][:, 0][mask_nans] + 2454833.0,
                   fluxes=hdf5_file[str(kic)][:, 1][mask_nans],
                   errors=hdf5_file[str(kic)][:, 2][mask_nans],
                   quarters=hdf5_file[str(kic)][:, 4][mask_nans],
                   name=name, params=params)

    def phases(self):
        params = self.params
        phase = ((self.times.jd - params.t0) % params.per)/params.per
        phase[phase > 0.5] -= 1.0
        return phase

    def plot(self, transit_params=None, ax=None, quarter=None, show=True,
             phase=False, plot_kwargs={'color':'b', 'marker':'o', 'lw':0}):
        """
        Plot light curve.

        Parameters
        ----------
        transit_params : `~batman.TransitParams` (optional)
            Transit light curve parameters. Required if `phase` is `True`.
        ax : `~matplotlib.axes.Axes` (optional)
            Axis to make plot on top of
        quarter : float (optional)
            Plot only this Kepler quarter
        show : bool
            If `True`, call `matplotlib.pyplot.show` after plot is made
        phase : bool
            If `True`, map times in JD to orbital phases, which requires
            that `transit_params` be input also.
        plot_kwargs : dict
            Keyword arguments to pass to `~matplotlib` calls.
        """
        if quarter is not None:
            if hasattr(quarter, '__len__'):
                mask = np.zeros_like(self.fluxes).astype(bool)
                for q in quarter:
                    mask |= self.quarters == q
            else:
                mask = self.quarters == quarter
        else:
            mask = np.ones_like(self.fluxes).astype(bool)

        if ax is None:
            ax = plt.gca()

        if phase:
            x = (self.times.jd - transit_params.t0)/transit_params.per % 1
            x[x > 0.5] -= 1
        else:
            x = self.times.jd

        ax.plot(x[mask], self.fluxes[mask],
                **plot_kwargs)
        ax.set(xlabel='Time' if not phase else 'Phase',
               ylabel='Flux', title=self.name)

        if show:
            plt.show()

    def normalize_each_quarter(self, rename=None, polynomial_order=2,
                               plots=False):
        """
        Use polynomial fit to each quarter to normalize the data.

        Parameters
        ----------
        rename : str (optional)
            New name of the light curve after normalization
        polynomial_order : int (optional)
            Order of polynomial to fit to the out-of-transit fluxes. Default
            is 2.
        plots : bool (optional)
            Show diagnostic plots after normalization.
        """
        quarter_inds = list(set(self.quarters))
        quarter_masks = [quarter == self.quarters for quarter in quarter_inds]

        for quarter_mask in quarter_masks:

            polynomial = np.polyfit(self.times[quarter_mask].jd,
                                    self.fluxes[quarter_mask], polynomial_order)
            scaling_term = np.polyval(polynomial, self.times[quarter_mask].jd)
            self.fluxes[quarter_mask] /= scaling_term
            self.errors[quarter_mask] /= scaling_term

            if plots:
                plt.plot(self.times[quarter_mask], self.fluxes[quarter_mask])
                plt.show()

        if rename is not None:
            self.name = rename

    def delete_outliers(self):

        d = np.diff(self.fluxes)
        spikey = np.abs(d - np.median(d)) > 2.5*np.std(d)
        neighboring_spikes = spikey[1:] & spikey[:-1]
        opposite_signs = np.sign(d[1:]) != np.sign(d[:-1])
        outliers = np.argwhere(neighboring_spikes & opposite_signs) + 1
        #print('number bad fluxes: {0}'.format(len(outliers)))

        self.times = Time(np.delete(self.times.jd, outliers), format='jd')
        self.fluxes = np.delete(self.fluxes, outliers)
        self.errors = np.delete(self.errors, outliers)
        self.quarters = np.delete(self.quarters, outliers)

    def mask_out_of_transit(self, oot_duration_fraction=0.25,
                            flip=False):
        """
        Mask out the out-of-transit light curve based on transit parameters

        Parameters
        ----------
        oot_duration_fraction : float (optional)
            Fluxes from what fraction of a transit duration of the
            out-of-transit light curve should be included in the mask?
        flip : bool (optional)
            If `True`, mask in-transit rather than out-of-transit.

        Returns
        -------
        d : dict
            Inputs for a new `LightCurve` object with the mask applied.
        """
        # Fraction of one duration to capture out of transit
        params = self.params
        phased = (self.times.jd - params.t0) % params.per
        near_transit = ((phased < params.duration*(0.5 + oot_duration_fraction)) |
                        (phased > params.per - params.duration*(0.5 + oot_duration_fraction)))
        if flip:
            near_transit = ~near_transit
        sort_by_time = np.argsort(self.times[near_transit].jd)
        return dict(times=self.times[near_transit][sort_by_time],
                    fluxes=self.fluxes[near_transit][sort_by_time],
                    errors=self.errors[near_transit][sort_by_time],
                    quarters=self.quarters[near_transit][sort_by_time],
                    params=self.params)

    def mask_in_transit(self, oot_duration_fraction=0.25):
        """
        Mask out the in-transit light curve based on transit parameters

        Parameters
        ----------
        oot_duration_fraction : float (optional)
            Fluxes from what fraction of a transit duration of the
            out-of-transit light curve should be included in the mask?

        Returns
        -------
        d : dict
            Inputs for a new `LightCurve` object with the mask applied.
        """
        params = self.params
        return self.mask_out_of_transit(params, flip=True,
                                        oot_duration_fraction=oot_duration_fraction)

    def get_transit_light_curves(self, plots=False):
        """
        For a light curve with transits only (i.e. like one returned by
        `LightCurve.mask_out_of_transit`), split up the transits into their
        own light curves, return a list of `TransitLightCurve` objects.

        Parameters
        ----------
        plots : bool
            Make diagnostic plots.

        Returns
        -------
        transit_light_curves : list
            List of `TransitLightCurve` objects
        """
        params = self.params
        time_diffs = np.diff(sorted(self.times.jd))
        diff_between_transits = params.per/2.
        split_inds = np.argwhere(time_diffs > diff_between_transits) + 1

        if len(split_inds) > 1:

            split_ind_pairs = [[0, split_inds[0][0]]]
            split_ind_pairs.extend([[split_inds[i][0], split_inds[i+1][0]]
                                     for i in range(len(split_inds)-1)])
            split_ind_pairs.extend([[split_inds[-1], len(self.times)]])

            transit_light_curves = []
            counter = -1
            for start_ind, end_ind in split_ind_pairs:
                counter += 1
                if plots:
                    plt.plot(self.times.jd[start_ind:end_ind],
                             self.fluxes[start_ind:end_ind], '.-')
                if hasattr(start_ind, '__len__'):
                    start_ind = start_ind[0]
                parameters = dict(times=self.times[start_ind:end_ind],
                                  fluxes=self.fluxes[start_ind:end_ind],
                                  errors=self.errors[start_ind:end_ind],
                                  quarters=self.quarters[start_ind:end_ind],
                                  name=counter,
                                  params=self.params)
                transit_light_curves.append(TransitLightCurve(**parameters))
            if plots:
                plt.show()
        else:
            transit_light_curves = []

        return transit_light_curves

    def get_available_quarters(self):
        """
        Get which quarters are available in this `LightCurve`

        Returns
        -------
        qs : list
            List of unique quarters available.
        """
        return list(set(self.quarters))

    def get_quarter(self, quarter):
        """
        Get a copy of the data from within `LightCurve` during one Kepler
        quarter.

        Parameters
        ----------
        quarter : int
            Kepler Quarter

        Returns
        -------
        lc : `LightCurve`
            Light curve from one Kepler Quarter
        """
        this_quarter = self.quarters == quarter
        return LightCurve(times=self.times[this_quarter],
                          fluxes=self.fluxes[this_quarter],
                          errors=self.errors[this_quarter],
                          quarters=self.quarters[this_quarter],
                          name=self.name + '_quarter_{0}'.format(quarter))

    @property
    def times_jd(self):
        """
        Get the times in this light curve in JD.

        Returns
        -------
        t_jd : `~numpy.ndarray`
            Julian dates.
        """
        return self.times.jd

    def split_at_index(self, index):
        """
        Split the light curve into two light curves, at ``index``
        """
        return (LightCurve(times=self.times[:index], fluxes=self.fluxes[:index],
                           errors=self.errors[:index], quarters=self.quarters[:index],
                           name=self.name),
                LightCurve(times=self.times[index:], fluxes=self.fluxes[index:],
                           errors=self.errors[index:], quarters=self.quarters[index:],
                           name=self.name))

    def transit_model(self, short_cadence=False):
        transit_params = self.params
        # (1 * u.min).to(u.day).value
        if short_cadence:
            exp_time = (1 * u.min).to(u.day).value #(6.019802903 * 10 * u.s).to(u.day).value
            supersample = 10
        else:
            exp_time = (6.019802903 * 10 * 30 * u.s).to(u.day).value
            supersample = 10

        m = batman.TransitModel(transit_params, self.times.jd,
                                supersample_factor=supersample,
                                exp_time=exp_time)
        model_flux = m.light_curve(transit_params)
        return model_flux


class TransitLightCurve(LightCurve):
    """
    Container for a single transit light curve. Subclass of `LightCurve`.
    """
    def __init__(self, times=None, fluxes=None, errors=None, quarters=None,
                 name=None, params=None):
        """
        Parameters
        ----------
        times : `~numpy.ndarray`
            Times in JD
        fluxes : `~numpy.ndarray`
            Fluxes (normalized or not)
        errors : `~numpy.ndarray`
            Uncertainties on the fluxes
        quarters : `~numpy.ndarray` (optional)
            Kepler Quarter for each flux
        name : str
            Name this light curve (optional)
        """

        if isinstance(times[0], Time) and isinstance(times, np.ndarray):
            times = Time(times)
        elif not isinstance(times, Time):
            times = Time(times, format='jd')
        self.times = times
        self.fluxes = fluxes
        self.errors = errors
        if self.times is not None and quarters is None:
            quarters = np.zeros_like(self.fluxes) - 1
        self.quarters = quarters
        self.name = name if type(name) is str else str(name)
        self.rescaled = False
        self.params = params

    def fit_linear_baseline(self, cadence=30*u.min,
                            return_near_transit=False, plots=False):
        """
        Find OOT portions of transit light curve using similar method to
        `LightCurve.mask_out_of_transit`, fit linear baseline to OOT.

        Parameters
        ----------
        cadence : `~astropy.units.Quantity` (optional)
            Length of the exposure time for each flux. Default is 1 min.
        return_near_transit : bool (optional)
            Return the mask for times in-transit.

        Returns
        -------
        linear_baseline : `numpy.ndarray`
            Baseline trend of out-of-transit fluxes
        near_transit : `numpy.ndarray` (optional)
            The mask for times in-transit.
        """
        params = self.params
        cadence_buffer = cadence.to(u.day).value
        get_oot_duration_fraction = 0
        phased = (self.times.jd - params.t0) % params.per
        near_transit = ((phased < params.duration *
                         (0.5 + get_oot_duration_fraction) + cadence_buffer) |
                        (phased > params.per - params.duration *
                         (0.5 + get_oot_duration_fraction) - cadence_buffer))

        # Remove linear baseline trend
        order = 1
        linear_baseline = np.polyfit(self.times.jd[~near_transit],
                                     self.fluxes[~near_transit], order)
        linear_baseline_fit = np.polyval(linear_baseline, self.times.jd)

        if plots:
            fig, ax = plt.subplots(1, 2, figsize=(15,6))
            ax[0].axhline(1, ls='--', color='k')
            ax[0].plot(self.times.jd, linear_baseline_fit, 'r')
            ax[0].plot(self.times.jd, self.fluxes, 'bo')
            plt.show()

        if return_near_transit:
            return linear_baseline, near_transit
        else:
            return linear_baseline

    def remove_linear_baseline(self, plots=False, cadence=30*u.min):
        """
        Find OOT portions of transit light curve using similar method to
        `LightCurve.mask_out_of_transit`, fit linear baseline to OOT,
        divide whole light curve by that fit.

        Parameters
        ----------
        cadence : `~astropy.units.Quantity` (optional)
            Length of the exposure time for each flux. Default is 1 min.
        plots : bool (optional)
            Show diagnostic plots.
        """
        params = self.params
        linear_baseline, near_transit = self.fit_linear_baseline(cadence=cadence,
                                                                 return_near_transit=True)
        linear_baseline_fit = np.polyval(linear_baseline, self.times.jd)
        self.fluxes =  self.fluxes/linear_baseline_fit
        self.errors = self.errors/linear_baseline_fit

        if plots:
            fig, ax = plt.subplots(1, 2, figsize=(15,6))
            ax[0].axhline(1, ls='--', color='k')
            ax[0].plot(self.times.jd, self.fluxes, 'o')
            ax[0].set_title('before trend removal')

            ax[1].set_title('after trend removal')
            ax[1].axhline(1, ls='--', color='k')
            ax[1].plot(self.times.jd, self.fluxes, 'o')
            plt.show()

    def remove_polynomial_baseline(self, order=2, plots=False, cadence=30*u.min):
        """
        Find OOT portions of transit light curve using similar method to
        `LightCurve.mask_out_of_transit`, fit linear baseline to OOT,
        divide whole light curve by that fit.

        Parameters
        ----------
        cadence : `~astropy.units.Quantity` (optional)
            Length of the exposure time for each flux. Default is 1 min.
        plots : bool (optional)
            Show diagnostic plots.
        """
        params = self.params
        poly_baseline, near_transit = self.fit_polynomial_baseline(order=order,
                                                                   cadence=cadence,
                                                                   return_near_transit=True)
        poly_baseline_fit = np.polyval(poly_baseline, self.times.jd)
        self.fluxes = self.fluxes / poly_baseline_fit
        self.errors = self.errors / poly_baseline_fit

        if plots:
            fig, ax = plt.subplots(1, 2, figsize=(15,6))
            ax[0].axhline(1, ls='--', color='k')
            ax[0].plot(self.times.jd, self.fluxes, 'o')
            ax[0].set_title('before trend removal')

            ax[1].set_title('after trend removal')
            ax[1].axhline(1, ls='--', color='k')
            ax[1].plot(self.times.jd, self.fluxes, 'o')
            plt.show()

    def scale_by_baseline(self, linear_baseline_params):
        if not self.rescaled:
            scaling_vector = np.polyval(linear_baseline_params, self.times.jd)
            self.fluxes *= scaling_vector
            self.errors *= scaling_vector
            self.rescaled = True

    def fit_polynomial_baseline(self, order=2, cadence=30*u.min,
                                plots=False, mask=None):
        """
        Find OOT portions of transit light curve using similar method to
        `LightCurve.mask_out_of_transit`, fit linear baseline to OOT
        """
        params = self.params
        if mask is None:
            mask = np.ones(len(self.fluxes)).astype(bool)
        cadence_buffer = cadence.to(u.day).value
        get_oot_duration_fraction = 0
        phased = (self.times.jd[mask] - params.t0) % params.per
        near_transit = ((phased < params.duration*(0.5 + get_oot_duration_fraction) + cadence_buffer) |
                        (phased > params.per - params.duration*(0.5 + get_oot_duration_fraction) - cadence_buffer))

        # Remove polynomial baseline trend after subtracting the times by its
        # mean -- this improves numerical stability for polyfit

        downscaled_times = self.times.jd - self.times.jd.mean()

        if len(downscaled_times[mask][~near_transit]) > 0:
            polynomial_baseline = np.polyfit(downscaled_times[mask][~near_transit],
                                             self.fluxes[mask][~near_transit],
                                             order)
            polynomial_baseline_fit = np.polyval(polynomial_baseline,
                                                 downscaled_times)
        else:
            polynomial_baseline_fit = np.ones_like(near_transit)

        if plots:
            fig, ax = plt.subplots(1, 2, figsize=(15,6))
            ax[0].axhline(1, ls='--', color='k')
            ax[0].plot(self.times.jd, polynomial_baseline_fit, 'r')
            ax[0].plot(self.times.jd, self.fluxes, 'bo')
            if mask is not None:
                ax[0].plot(self.times.jd[~mask], self.fluxes[~mask], 'ro')
            plt.show()

        return polynomial_baseline_fit

    def subtract_polynomial_baseline(self, plots=False, order=2,
                                     cadence=30*u.min):
        """
        Find OOT portions of transit light curve using similar method to
        `LightCurve.mask_out_of_transit`, fit polynomial baseline to OOT,
        subtract whole light curve by that fit.
        """
        params = self.params
        polynomial_baseline_fit = self.fit_polynomial_baseline(cadence=cadence,
                                                               order=order)
        self.fluxes = self.fluxes - polynomial_baseline_fit
        self.errors = self.errors

        if plots:
            fig, ax = plt.subplots(1, 2, figsize=(15,6))
            ax[0].axhline(1, ls='--', color='k')
            ax[0].plot(self.times.jd, self.fluxes, 'o')
            #ax[0].plot(self.times.jd[near_transit], self.fluxes[near_transit], 'ro')
            ax[0].set_title('before trend removal')

            ax[1].set_title('after trend removal')
            ax[1].axhline(1, ls='--', color='k')
            ax[1].plot(self.times.jd, self.fluxes, 'o')
            plt.show()


    def subtract_add_divide_without_outliers(self, quarterly_max,
                                             order=2, cadence=30*u.min,
                                             outlier_error_multiplier=50,
                                             outlier_tolerance_depth_factor=0.20,
                                             plots=False):

        if len(self.times_jd) > 0:
            init_baseline_fit = self.fit_polynomial_baseline(order=order,
                                                             cadence=cadence)

            # Subtract out a transit model
            transit_model = generate_lc_depth(self.times_jd, self.params.rp**2, self.params)

            lower_outliers = (transit_model*init_baseline_fit - self.fluxes >
                              self.fluxes.mean() * outlier_tolerance_depth_factor *
                              self.params.rp**2)

            self.errors[lower_outliers] *= outlier_error_multiplier

            final_baseline_fit = self.fit_polynomial_baseline(order=order,
                                                              cadence=cadence,
                                                              mask=~lower_outliers)

            self.fluxes = self.fluxes - final_baseline_fit
            self.fluxes += quarterly_max
            self.fluxes /= quarterly_max
            self.errors /= quarterly_max

            if plots:
                plt.errorbar(self.times.jd, self.fluxes, self.errors, fmt='o')
                plt.plot(self.times.jd[lower_outliers],
                         self.fluxes[lower_outliers], 'rx')
                plt.show()

    def chi2_lc_3param(self, p):
        """
        Compute chi^2 for a three-parameter light curve model with parameters
        in tuple ``p`` consisting of the planet-star radius ratio, and two
        limb-darkening parameters (Kipping 2013)
        """
        rp, q1, q2 = p
        model = lc_3param(self.times_jd, rp, q2, q2, self.params)

        mask_nans = np.logical_not(np.isnan(self.fluxes) |
                                   np.isnan(self.errors))

        return np.sum((model[mask_nans] - self.fluxes[mask_nans])**2 /
                      (2*self.errors[mask_nans])**2)

    def fit_lc_3param(self):
        """
        Fit three-parameter light curve model, and replace the transit parameters
        in ``self.params`` with the best fit planet-star radius ratio, and two
        limb-darkening parameters (Kipping 2013).
        """
        q1_init, q2_init = u2q(*self.params.u)
        initp = [self.params.rp, q1_init, q2_init]

        bounds = [(0.1 * self.params.rp, 2 * self.params.rp), (0, 1), (0, 1)]
        results = fmin_l_bfgs_b(self.chi2_lc_3param, initp, approx_grad=True,
                                bounds=bounds, iprint=0)
        bestp = results[0]
        best_rp, best_q1, best_q2 = bestp
        self.params.rp = best_rp
        self.params.u = q2u(best_q1, best_q2)


def concatenate_transit_light_curves(light_curve_list, name=None):
    """
    Combine multiple transit light curves into one `TransitLightCurve` object.

    Parameters
    ----------
    light_curve_list : list
        List of `TransitLightCurve` objects
    name : str
        Name of new light curve

    Returns
    -------
    tlc : `TransitLightCurve`
        Concatenated transit light curves
    """
    times = []
    fluxes = []
    errors = []
    quarters = []
    for light_curve in light_curve_list:
        times.append(light_curve.times.jd)
        fluxes.append(light_curve.fluxes)
        errors.append(light_curve.errors)
        quarters.append(light_curve.quarters)
    times, fluxes, errors, quarters = [np.concatenate(i)
                                       for i in [times, fluxes,
                                                 errors, quarters]]

    times = Time(times, format='jd')
    return TransitLightCurve(times=times, fluxes=fluxes, errors=errors,
                             quarters=quarters, name=name, params=light_curve.params)


def concatenate_light_curves(light_curve_list, name=None):
    """
    Combine multiple transit light curves into one `TransitLightCurve` object.

    Parameters
    ----------
    light_curve_list : list
        List of `TransitLightCurve` objects
    name : str
        Name of new light curve

    Returns
    -------
    tlc : `TransitLightCurve`
        Concatenated transit light curves
    """
    times = []
    fluxes = []
    errors = []
    quarters = []
    for light_curve in light_curve_list:
        times.append(light_curve.times.jd)
        fluxes.append(light_curve.fluxes)
        errors.append(light_curve.errors)
        quarters.append(light_curve.quarters)
    times, fluxes, errors, quarters = [np.concatenate(i)
                                       for i in [times, fluxes,
                                                 errors, quarters]]

    times = Time(times, format='jd')
    return LightCurve(times=times, fluxes=fluxes, errors=errors,
                      quarters=quarters, name=name, params=light_curve.params[0])


def subtract_add_divide(whole_lc, transits):
    """
    Normalize transit light curves in ``transits`` with the
    "subtract-add-divide" method.

    Parameters
    ----------
    whole_lc : `~salter.LightCurve` or subclass
        Light curve over the whole kepler mission
    transits : list of `~salter.LightCurve` or subclasses
        Transit light curves
    """
    # Compute maxes for each quarter
    available_quarters = whole_lc.get_available_quarters()
    quarters = [whole_lc.get_quarter(q) for q in whole_lc.get_available_quarters()]

    quarterly_maxes = {}
    set_upper_limit = 4e10
    for i, quarter_number, lc in zip(range(len(available_quarters)),
                                     available_quarters, quarters):
        fluxes = lc.fluxes[lc.fluxes < set_upper_limit]
        smoothed_fluxes = gaussian_filter(fluxes, sigma=20)
        quarterly_maxes[quarter_number] = np.max(smoothed_fluxes)

    for transit in transits:
        transit.subtract_add_divide_without_outliers(quarterly_max=quarterly_maxes[lc.quarters[0]],
                                                     plots=False)
