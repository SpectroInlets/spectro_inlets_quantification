# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains functionality for processing of raw mass-spec i.e. peak fitting and
background subtraction

Variables with abbreviated or non-descriptive names, e.g. physical quantities:
* S or y: signal in [A]
* M or x: m/z ratio in atomic units [atomic mass unit / fundamental charge]
* bg: background in [A]
* t: time in [s]

"""
from typing import Union, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit

from .constants import STANDARD_WIDTH
from .tools import make_axis
from .exceptions import PeakFitError


# ---------- Generic Peak -------------- #

BackGround = Union[NDArray, float]


class Peak:
    def __init__(
        self,
        x: NDArray = None,
        y: NDArray = None,
        t: float = None,
        bg: BackGround = None,
        error: bool = False,
    ):
        """Initiate the peak with x, y and bg

        Args:
            x (np.array): m/z values in spectrum
            y (np.array): Signal in spectrum in [A]. Vector with same length as x.
            t (float): The measurement time of the peak in [s]
            bg (np.array or float): signal corresponding to background in [A]. Can be a
                vector with the same length as x and y, or a constant. The latter makes
                use of broadcasting to subtract bg from y.
                FIXME: bg should probably be removed from Peak (now it's just 0), since
                    signal.SignalProcessor takes care of background (also 0 for now).
            error (bool): whether the Peak is known to have or result from an error
                This by default means that it has zero signal.
        """
        self.x = x
        self.y_raw = y  # y before background subtraction
        if bg is not None:
            self.y = y - bg
            self.bg = bg
        else:
            self.y = y
            self.bg = None
        self.t = t
        self._height: Optional[float] = None
        self._center: Optional[float] = None
        self._width: Optional[float] = None
        self._y_fit = None
        self._x_fit = None
        self._integral = None
        self.error = error  # used by PeakSeries and SignalProcessor

    def reset_bg(self, bg: BackGround) -> None:
        """Reset the background, calculate fresh background-subtracted y"""
        self.bg = bg
        self.y = self.y_raw - bg

    def background_function(self, x):
        """Return the background as float (vector) for the m/z value (values) in x"""
        bg = self.bg if self.bg is not None else 0
        try:
            iter(bg)
        except TypeError:
            try:
                shape = x.shape
            except AttributeError:
                return bg
            else:
                return np.ones(shape) * bg
        else:
            return np.interp(x, self.x, bg)

    @property
    def height(self) -> float:
        """The height of the peak in [A]"""
        if not self._height:
            self._height = self.calc_height()
        return self._height

    @property
    def center(self) -> float:
        """The center of the peak in m/z units"""
        if not self._center:
            self._center = self.calc_center()
        return self._center

    @property
    def width(self) -> float:
        """The width of the peak in m/z units"""
        if not self._width:
            self._width = self.calc_width()
        return self._width

    @property
    def integral(self):
        if not self._integral:
            self.calc_integral()
        return self._integral

    @property
    def y_fit(self):
        return self._y_fit

    @property
    def x_fit(self):
        """dummy x-vector only used for plotting the fit. 40 ppamu resolution."""
        if self._x_fit is None:
            self._x_fit = np.arange(self.x[0], self.x[-1], 0.025)
        return self._x_fit

    def calc_height(self) -> float:
        """Simplest implementation of height, return maximum signal value in [A].

        This should be overwritten in inheriting classes like GaussPeak
        """
        return max(self.y)

    def calc_center(self) -> float:
        """Simplest implementation of center, return m/z at maximum signal value.

        This should be overwritten in inheriting classes like GaussPeak"""
        return self.x[np.argmax(self.y)]

    def calc_width(self) -> float:
        """Return diff. between m/z of first and last point where signal > 10% of height

        This should be overwritten in inheriting classes like GaussPeak"""
        x_in_peak = self.x[self.y > self.height / 10]
        return x_in_peak[-1] - x_in_peak[0]

    def calc_integral(self):
        """Return trapazoidal integral of Peak in [m/z * A]"""
        return np.trapz(self.y, self.x)

    def calc_signal(self):
        """The signal is by default the height of the peak in [A]."""
        if self.error:
            return 0
        return self.height

    @property
    def signal(self):
        """The signal is by default the height of the peak in [A]."""
        return self.calc_signal()

    def plot(self, ax="new", **kwargs):
        """Plot the peak and its background

        dev-time only

        Args:
            ax (Axis or str): axis to plot on. "new" makes new axis
            kwargs: additional key-word args are taken as plotting spec
        """
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("m/z")
            ax.set_ylabel("signal / [A]")
        ax.plot(self.x, self.y_raw, **kwargs)
        if self.bg is not None:
            bg_plot = self.background_function(self.x)
            bg_kwargs = {**kwargs, "linestyle": "--"}
            ax.plot(self.x, bg_plot, **bg_kwargs)
        if self.y_fit is not None:
            fit_kwargs = {**kwargs, "linestyle": ":", "label": "fit"}
            ax.plot(self.x_fit, self.y_fit, **fit_kwargs)
        return ax


# ---------- Gaussian peak -------------- #


def gauss_fun(x, center, sigma, height):
    """Return the Gaussian function of x with the given parameters"""
    y = height * np.exp(-((x - center) ** 2) / (2 * sigma**2))
    return y


class GaussPeak(Peak):
    """A Peak with a gauss fit"""

    def __init__(self, center=None, sigma=None, tolerance=0.2, *args, **kwargs):
        """When a GaussPeak is initiated, it fits x to y with a Gauss.

        Args:
             center (float): the center of the peak, if known, in m/z units
             sigma (float): the sigma of the peak, if known, in m/z units
             tolerance (float): the relative_square error above which a
                 PeakFitError is raised
             additional args and kwargs go to Peak.__init__()
        """
        super().__init__(*args, **kwargs)
        self._fwhm = None
        self._sigma = None
        self.tolerance = tolerance
        self.fit_gauss(center=center, sigma=sigma)

    def fit_gauss(self, center=None, sigma=None):
        """Find the Gauss function of self.x minimizing the square error wrt self.y

        Either or both of center and sigma can be given if they are known.
        If center and/or sigma are not known, they are fit.
        The height is always fit.
        The function sets the properties center, sigma=width, and height of the peak,
        as well as its integral (the analytical integral of the gauss) and y_fit.

        Args:
             center (float): the center of the peak, if known, in m/z units
             sigma (float): the sigma of the peak, if known, in m/z units

        Returns tuple: (center, sigma, height) of the fit peak
        """
        x, y = self.x, self.y

        guess_c = (x[-1] + x[0]) / 2
        guess_s = STANDARD_WIDTH / 3  # 1/3 converts width at 10% max to sigma
        # ^ TODO: replace this emperical guess with the analytical number
        guess_h = max(y)

        if center is not None and sigma is not None:
            # Then we just need to fit the height
            def gauss_i(x_, height_):
                return gauss_fun(x_, center=center, sigma=sigma, height=height_)

            guess = guess_h
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            height = popt[0]

        elif center is not None:

            def gauss_i(x_, sigma_, height_):
                return gauss_fun(x_, center=center, sigma=sigma_, height=height_)

            guess = [guess_s, guess_h]
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            sigma, height = popt[0], popt[1]
        elif sigma is not None:

            def gauss_i(x_, center_, height_):
                return gauss_fun(x_, center=center_, sigma=sigma, height=height_)

            guess = [guess_c, guess_h]
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            center, height = popt[0], popt[1]
        else:

            def gauss_i(x_, center_, sigma_, height_):
                return gauss_fun(x_, center=center_, sigma=sigma_, height=height_)

            guess = [guess_c, guess_s, guess_h]
            try:
                popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
                center, sigma, height = popt[0], popt[1], popt[2]
            except RuntimeError as e:
                raise PeakFitError(f"curve_fit raises RunTimeError('{e}')")

        sigma = abs(sigma)

        # print(f'center={center}, sigma={sigma}, height={height}') # debugging
        y_fit = gauss_fun(self.x_fit, center, sigma, height) + self.background_function(
            self.x_fit
        )
        integral_f = np.sqrt(2 * np.pi) * height * sigma
        self._center = center
        self._sigma = sigma  # note that GaussPeak defines its width as its sigma
        self._height = height

        rse = self.relative_square_error
        if rse > self.tolerance:
            raise PeakFitError(
                f"relative_square_error = {rse} > tolerance = {self.tolerance}"
            )

        self._y_fit = y_fit
        self._integral = integral_f  # ... and its integral as the analytical integral
        return center, sigma, height

    def y_of_x(self, x):
        """Return the fit background-subtracted signal given m/z"""
        y = gauss_fun(x, center=self.center, sigma=self.sigma, height=self.height)
        return y

    def y_raw_of_x(self, x):
        """Return the fit raw signal given m/z"""
        y = self.y_of_x(x)
        y_raw = y + self.background_function(x)
        return y_raw

    @property
    def relative_square_error(self):
        """Return the relative square error of the fit as a float"""
        error = self.y_of_x(self.x) - self.y
        return error.dot(error) / self.y.dot(self.y)

    @property
    def sigma(self):
        """The width for a GaussPeak is defined as its sigma"""
        return self._sigma

    @property
    def width(self):
        """The width for a GaussPeak is defined as its sigma"""
        return self._sigma

    @property
    def fwhm(self):
        """Analytical full width half maximum of the gauss peak"""
        fwhm = 2 * np.sqrt(2 * np.log(2)) * self.sigma
        return fwhm


PEAK_CLASSES = {"simple": Peak, "gauss": GaussPeak}
