# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains functionality for processing of raw mass-spec i.e. peak fitting and
background subtraction

Variables with abbreviated or non-descriptive names, e.g. physical quantities:
 * S or y: signal in [A]
 * M or x: m/z ratio in atomic units [atomic mass unit / fundamental charge]
 * bg: background in [A]
 * t: time in [s]

"""
from typing import Union, Optional, cast, TYPE_CHECKING, Any, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit  # type: ignore

from .constants import STANDARD_WIDTH
from .tools import make_axis
from .exceptions import PeakFitError
from .compatability import TypeAlias


if TYPE_CHECKING:
    from matplotlib import pyplot

# ---------- Generic Peak -------------- #


FLOAT_OR_FLOAT_ARRAY: TypeAlias = Union[NDArray[np.float64], float]
BACKGROUND: TypeAlias = FLOAT_OR_FLOAT_ARRAY
FLOAT_ARRAY: TypeAlias = NDArray[np.float64]


class Peak:
    """Class which implements a basic peak"""

    def __init__(
        self,
        x: FLOAT_ARRAY,
        y: FLOAT_ARRAY,
        t: Optional[float] = None,
        bg: Optional[BACKGROUND] = None,
        error: bool = False,
    ):
        """Initiate the peak with x, y and bg

        Args:
            x (np.array): m/z values in spectrum
            y (np.array): Signal of spectrum in [A]. Vector with same length as x.
            t (float): The measurement time of the peak in [s]
            bg (np.array or float): Signal corresponding to background in [A]. Can be a
                vector with the same length as x and y, or a constant. The latter makes
                use of broadcasting to subtract bg from y.
            error (bool): Whether the Peak is known to have or result from an error
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
        self._y_fit: Optional[FLOAT_ARRAY] = None
        self._x_fit: Optional[FLOAT_ARRAY] = None
        self._integral: Optional[float] = None
        self.error = error  # used by PeakSeries and SignalProcessor

    def background_function(self, x: NDArray) -> NDArray:
        """Return the background as float (vector) for the m/z value (values) in x"""
        bg = self.bg if self.bg is not None else 0.0
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
    def integral(self) -> float:
        """The integral of the signal in [m/z * A]"""
        if not self._integral:
            self._integral = self.calc_integral()
        return self._integral

    @property
    def y_fit(self) -> Optional[FLOAT_ARRAY]:
        """Returns an array with fitted y-values"""
        return self._y_fit

    @property
    def x_fit(self) -> FLOAT_ARRAY:
        """dummy x-vector only used for plotting the fit. 40 ppamu resolution."""
        if self._x_fit is None:
            self._x_fit = np.arange(self.x[0], self.x[-1], 0.025)
        return self._x_fit

    def calc_height(self) -> float:
        """Simplest implementation of height, return maximum signal value in [A].

        This should be overwritten in inheriting classes like GaussPeak
        """
        return cast(float, max(self.y))

    def calc_center(self) -> float:
        """Simplest implementation of center, return m/z at maximum signal value.

        This should be overwritten in inheriting classes like GaussPeak
        """
        return cast(float, self.x[np.argmax(self.y)])

    def calc_width(self) -> float:
        """Return diff. between m/z of first and last point where signal > 10% of height

        This should be overwritten in inheriting classes like GaussPeak"""
        x_in_peak = self.x[self.y > self.height / 10]
        return cast(float, x_in_peak[-1] - x_in_peak[0])

    def calc_integral(self) -> float:
        """Return trapezoidal integral of Peak in [m/z * A]"""
        return cast(float, np.trapz(self.y, self.x))

    def calc_signal(self) -> float:
        """The signal is by default the height of the peak in [A]"""
        if self.error:
            return 0
        return self.height

    @property
    def signal(self) -> float:
        """The signal is by default the height of the peak in [A]"""
        return self.calc_signal()

    def plot(
        self, ax: Union[str, "pyplot.Axes"] = "new", **kwargs: Any
    ) -> "pyplot.Axes":  # pragma: no cover
        """Plot the peak and its background

        dev-time only

        Args:
            ax (Axis or str): axis to plot on. "new" makes new axis
            kwargs: additional key-word args are taken as plotting spec

        Returns:
            Axis: The axis that was used for plotting
        """
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("m/z")
            ax.set_ylabel("signal / [A]")
        ax = cast("pyplot.Axes", ax)
        ax.plot(self.x, self.y_raw, **kwargs)
        if self.bg is not None:
            bg_plot = self.background_function(self.x)
            bg_kwargs = {**kwargs, "linestyle": "--"}
            ax.plot(self.x, bg_plot, **bg_kwargs)  # type: ignore
        if self.y_fit is not None:
            fit_kwargs = {**kwargs, "linestyle": ":", "label": "fit"}
            ax.plot(self.x_fit, self.y_fit, **fit_kwargs)  # type: ignore
        return ax


# ---------- Gaussian peak -------------- #


def gauss_fun(
    x: FLOAT_OR_FLOAT_ARRAY, center: float, sigma: float, height: float
) -> FLOAT_OR_FLOAT_ARRAY:
    """Return the Gaussian function of x with the given parameters"""
    y = height * np.exp(-((x - center) ** 2) / (2 * sigma**2))
    return y


class GaussPeak(Peak):
    """A Peak with a gauss fit"""

    def __init__(
        self,
        x: FLOAT_ARRAY,
        y: FLOAT_ARRAY,
        center: Optional[float] = None,
        sigma: Optional[float] = None,
        tolerance: Optional[float] = 0.2,
        **kwargs: Any,
    ) -> None:
        """When a GaussPeak is initiated, it fits x to y with a Gauss.

        Args:
            x (np.array): m/z values in spectrum
            y (np.array): Signal of spectrum in [A]. Vector with same length as x.
            center (float): The center of the peak, if known, in m/z units
            sigma (float): The sigma of the peak, if known, in m/z units
            tolerance (float): The relative_square error above which a
                PeakFitError is raised
            kwargs (Any): Additional kwargs go to `Peak.__init__`, most notably the t and bg
                arguments
        """
        super().__init__(x, y, **kwargs)
        self._fwhm: Optional[float] = None
        self.sigma: Optional[float] = None
        self.tolerance = tolerance
        self.fit_gauss(center=center, sigma=sigma)

    def fit_gauss(
        self, center: Optional[float] = None, sigma: Optional[float] = None
    ) -> Tuple[float, float, float]:
        """Find the Gauss function of self.x minimizing the square error with regards to self.y

        Either or both of center and sigma can be given if they are known. If center and/or
        sigma are not known, they are fitted. The height is always fitted. The function sets
        the properties center, sigma=width, and height of the peak, as well as its integral
        (the analytical integral of the gauss) and y_fit.

        Args:
             center (float): The center of the peak, if known, in m/z units
             sigma (float): The sigma of the peak, if known, in m/z units

        Returns tuple: (center, sigma, height) of the fit peak
        """
        x, y = self.x, self.y

        guess_c = (x[-1] + x[0]) / 2
        guess_s = STANDARD_WIDTH / 3  # 1/3 converts width at 10% max to sigma
        # ^ TODO: replace this empirical guess with the analytical number
        guess_h = max(y)

        if center is not None and sigma is not None:
            # Then we just need to fit the height
            def gauss_i_height_only(x_: float, height_: float) -> FLOAT_OR_FLOAT_ARRAY:
                return gauss_fun(x_, center=center, sigma=sigma, height=height_)

            guess = guess_h
            popt, pcov = curve_fit(gauss_i_height_only, x, y, p0=guess)
            height = popt[0]

        elif center is not None:

            def gauss_i_sigma_height(
                x_: float, sigma_: float, height_: float
            ) -> FLOAT_OR_FLOAT_ARRAY:
                return gauss_fun(x_, center=center, sigma=sigma_, height=height_)

            guess = [guess_s, guess_h]
            popt, pcov = curve_fit(gauss_i_sigma_height, x, y, p0=guess)
            sigma, height = popt[0], popt[1]
        elif sigma is not None:

            def gauss_i_center_height(
                x_: float, center_: float, height_: float
            ) -> FLOAT_OR_FLOAT_ARRAY:
                return gauss_fun(x_, center=center_, sigma=sigma, height=height_)

            guess = [guess_c, guess_h]
            popt, pcov = curve_fit(gauss_i_center_height, x, y, p0=guess)
            center, height = popt[0], popt[1]
        else:

            def gauss_i_all(
                x_: float, center_: float, sigma_: float, height_: float
            ) -> FLOAT_OR_FLOAT_ARRAY:
                return gauss_fun(x_, center=center_, sigma=sigma_, height=height_)

            guess = [guess_c, guess_s, guess_h]
            try:
                popt, pcov = curve_fit(gauss_i_all, x, y, p0=guess)
                center, sigma, height = popt[0], popt[1], popt[2]
            except RuntimeError as e:
                raise PeakFitError(f"curve_fit raises RunTimeError('{e}')")

        sigma = abs(sigma)

        y_fit = gauss_fun(self.x_fit, center, sigma, height) + self.background_function(
            self.x_fit
        )
        integral_f = np.sqrt(2 * np.pi) * height * sigma
        self._center = center
        self.sigma = sigma  # note that GaussPeak defines its width as its sigma
        self._height = height

        rse = self.relative_square_error
        if rse > self.tolerance:
            raise PeakFitError(
                f"relative_square_error = {rse} > tolerance = {self.tolerance}"
            )

        self._y_fit = y_fit
        self._integral = integral_f  # ... and its integral as the analytical integral
        return center, sigma, height

    def y_of_x(self, x: FLOAT_ARRAY) -> FLOAT_ARRAY:
        """Return the fit background-subtracted signal given m/z"""
        y = gauss_fun(x, center=self.center, sigma=self.sigma, height=self.height)
        return y

    def y_raw_of_x(self, x: FLOAT_ARRAY) -> FLOAT_ARRAY:
        """Return the fit raw signal given m/z"""
        y = self.y_of_x(x)
        y_raw = y + self.background_function(x)
        return y_raw

    @property
    def relative_square_error(self) -> float:
        """Return the relative square error of the fit as a float"""
        error = self.y_of_x(self.x) - self.y
        return cast(float, error.dot(error) / self.y.dot(self.y))

    @property
    def width(self) -> float:
        """Return the width of the GaussPeak, defined as its sigma"""
        return self.sigma

    @property
    def fwhm(self) -> float:
        """Return the analytical full width half maximum of the gauss peak"""
        fwhm = 2 * np.sqrt(2 * np.log(2)) * self.sigma
        return cast(float, fwhm)


PEAK_CLASSES = {"simple": Peak, "gauss": GaussPeak}
