# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Converting raw MS data to one S_vec per mass, corrected for bg and non-linearity

The signal module is rather high up in the hierarchy so that it can be aware of
system conditions (Medium) and even the calibration (which I think is the right place
to store nonlinear coefficients).

The SignalProcessor, with help from the other classes of this and the peak module does:

* On initiation:

  * Load one or more background spectra
  * Load one or more linearity correction functions
* During operation:

  * Process channel data to produce signal values (S_M) to feed to quant.calibration. This includes:

    - Subtracting background from the channel data
    - Fitting the resulting peaks
  * Maybe adjusts backgrounds on the fly to match the peak-free areas of the measured full spectra?

The `SignalDict` and ``PeakSeries`` classes of this module are the (advanced) MID and spectrum
data classes of quant. They have useful data visualization methods.

PeakSeries is the main analysis tool of recorded datasets (a parser should load raw data
to a PeakSeries as the entry point to quant). It owns a SignalProcessor and uses it
to generate the data in its SignalDict.

Likewise, a live data analyzer can use SignalProcessor to fill a SignalDict continuously
during a measurement.
"""

from pathlib import Path  # noqa
import json
from typing import (
    Optional,
    Dict,
    Iterable,
    Tuple,
    Sequence,
    Generator,
    Union,
    cast,
    TYPE_CHECKING,
    Any,
    Type,
)

import numpy as np
from numpy.typing import NDArray
import time
from .constants import STANDARD_COLORS
from .custom_types import MASS_TO_SIGNAL
from .tools import mass_to_M, mass_to_pure_mass, make_axis
from .exceptions import PeakFitError
from .medium import Medium
from .peak import PEAK_CLASSES, Peak
from .config import Config
from .compatability import TypeAlias

if TYPE_CHECKING:
    from matplotlib import pyplot

CONFIG = Config()

FLOAT_ARRAY: TypeAlias = NDArray[np.float64]
MASS_LIST: TypeAlias = Iterable[str]
MASS_TO_ARRAY: TypeAlias = Dict[str, FLOAT_ARRAY]
TIME_SIGNAL_ARRAY_PAIR: TypeAlias = Tuple[FLOAT_ARRAY, FLOAT_ARRAY]
MASS_TO_TIME_SIGNAL_ARRAY_PAIR: TypeAlias = Dict[str, TIME_SIGNAL_ARRAY_PAIR]
TSPAN: TypeAlias = Sequence[float]
MSPAN: TypeAlias = Tuple[float, float]


class SignalDict:
    """Signals by (advanced) MID with memory. Index with mass for latest signal in [A].

    SignalDict is a wrapper around a set of mass-indexed signal vectors (S_vecs) and the
    corresponding time vectors (t_vecs) specifying when each signal was recorded. It
    is designed to be dynamic, making it easy to add MID signals as they are measured
    or calculated.

    Whenever a SignalDict is given a value assignment, it records the signal by
    appending it to ``S_vecs[mass]`` and appending the time to ``t_vecs[mass]``, removing
    the earliest one if needed due to the memory limitation set by max_depth.
    """

    def __init__(
        self,
        tstamp: Optional[float] = None,
        max_depth: Optional[int] = None,
        verbose: bool = False,
    ):
        """Initiate a SignalDict, optionally with initial signals

        Args:
            tstamp (float): The unix time defining the start of the measurement
            max_depth (int or None): The maximum length of the vectors, to save memory
            signals (dict): {mass: S_M_0} initial signals as float in [A] for each mass
            verbose (bool): Whether to print stuff to the terminal
        """
        self.tstamp = tstamp or time.time()
        self.S_vecs: MASS_TO_ARRAY = {}
        self.t_vecs: MASS_TO_ARRAY = {}
        self.max_depth = max_depth
        self.verbose = verbose

    @property
    def signals(self) -> MASS_TO_SIGNAL:
        """This is the dictionary with the latest signal for each mass"""
        S_dict = {mass: S_vec[-1] for mass, S_vec in self.S_vecs.items()}
        return S_dict

    def set_signal(self, mass: str, S_M: float, t: Optional[float] = None) -> None:
        """Add the signal S_M at mass with measurement t (by default, now)"""
        t = t if t is not None else (time.time() - self.tstamp)
        so_far = self.S_vecs.get(mass, np.array([]))
        t_so_far = self.t_vecs.get(mass, np.array([]))
        S_vec = np.append(so_far, S_M)
        t_vec = np.append(t_so_far, t)
        if self.max_depth and len(S_vec) > self.max_depth:
            S_vec = S_vec[-self.max_depth :]
            t_vec = t_vec[-self.max_depth :]
        self.S_vecs[mass] = S_vec
        self.t_vecs[mass] = t_vec

    def set_signals(self, signals: MASS_TO_SIGNAL, t: Optional[float] = None) -> None:
        """Set multiple signals at once with the same measurement time which defaults to now"""
        for mass, S_M in signals.items():
            self.set_signal(mass, S_M, t=t)

    def get_signal(
        self, mass: str, tspan: Optional[TSPAN] = None
    ) -> TIME_SIGNAL_ARRAY_PAIR:
        """Return (t, S_M), the time and signal vectors for mass during tspan"""
        t, S_M = self.t_vecs[mass], self.S_vecs[mass]
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t, S_M = t[mask], S_M[mask]
        return t, S_M

    def get_signals(
        self, mass_list: Optional[MASS_LIST] = None, tspan: Optional[TSPAN] = None
    ) -> MASS_TO_TIME_SIGNAL_ARRAY_PAIR:
        """Return the tspan-cut time and signal vecs for each mass in mass_list"""
        mass_list = mass_list or self.mass_list
        return {mass: self.get_signal(mass, tspan) for mass in mass_list}

    @property
    def mass_list(self) -> MASS_LIST:
        """The mass_list of a SignalDict lists the masses for which it has signals"""
        return list(self.S_vecs.keys())

    def signals_at(self, n: int) -> MASS_TO_SIGNAL:
        """Return the {mass: S_M} at the n'th iteration (negative is WRT most recent)"""
        return {mass: S_vec[n] for mass, S_vec in self.S_vecs.items()}

    def get_average_of_last(
        self,
        N: int = 50,
        mass_list: Optional[MASS_LIST] = None,
        t: Optional[float] = None,
    ) -> MASS_TO_SIGNAL:
        """Return signal averaged over last N iterations or t seconds for mass_list
        (default all)"""
        mass_list = mass_list or self.mass_list
        if t is not None:
            now = time.time() - self.tstamp
            tspan = [now - t, now]
            S_avg = {}
            for mass in mass_list:
                t_M, S_M = self.get_signal(mass, tspan)
                S_avg[mass] = np.mean(S_M)
        else:
            S_avg = {mass: np.mean(self.S_vecs[mass][-N:]) for mass in mass_list}
        return S_avg

    def __iter__(self) -> Generator[str, None, None]:
        """Iterating a SignalDict iterates over the masses for which it has signals"""
        for mass in self.mass_list:
            yield mass

    def __contains__(self, mass: str) -> bool:
        """A mass is in a SignalDict if it contains a signal at that mass"""
        return mass in self.mass_list

    def items(self) -> Generator[Tuple[str, FLOAT_ARRAY], None, None]:
        """An iterator over (mass, most recent advanced MID signal at mass)"""
        for mass in self:
            yield mass, self[mass]

    def __repr__(self) -> str:
        return f"SignalDict({self.signals})"

    def __getitem__(self, key: Union[str, int]) -> Union[float, FLOAT_ARRAY]:
        """Index with a mass to get its most recent signal.

        Index with a negative integer to get the signals at a past iteration.
        Use cinf-style indexing for vectors: "{mass}-y" for signal and "{mass}-x" for
        the corresponding time vector.
        """
        if cast(str, key) in self:
            return self.signals[cast(str, key)]
        elif isinstance(key, int):
            return self.signals_at(key)
        elif key.endswith("-x"):
            return self.t_vecs[key[:-2]]
        elif key.endswith("-y"):
            return self.S_vecs[key[:-2]]
        raise KeyError

    def __setitem__(self, mass: str, value: float) -> None:
        """Set the signal for mass to value with now as its measurement time"""
        self.set_signal(mass, value)

    def clear(self, mass: str) -> None:
        """Remove all knowledge of a mass from the SignalDict"""
        if mass in self:
            del self.S_vecs[mass]
            del self.t_vecs[mass]

    def clear_all(self) -> None:
        """Remove all knowledge at all masses from the SignalDict"""
        for mass in self:
            self.clear(mass)

    def plot(
        self,
        mass_list: MASS_LIST = None,
        tspan: TSPAN = None,
        ax: Union[str, "pyplot.Axes"] = "new",
    ) -> "pyplot.Axes":  # pragma: no cover
        """Plot the signal for each mass in mass_list, optionally restricted to tspan"""
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("time / [s]")
            ax.set_ylabel("signal / [A]")
        ax = cast("pyplot.Axes", ax)
        mass_list = mass_list or self.mass_list
        for mass in mass_list:
            color = STANDARD_COLORS.get(mass_to_pure_mass(mass), None)
            t, S = self.get_signal(mass, tspan)
            ax.plot(t, S, color=color, label=mass)
        return ax


class SignalProcessor:
    """Class for calculating one signal (S_M) from the MS data for each peak (y vs x)

    The x-and-y details is actually handled by the Peak class and its family via make_peak
    Indexing a SignalProcessor with the mass M gives S_M.
    """

    def __init__(
        self,
        *,
        mass_list: MASS_LIST = None,
        peak_type: str = "gauss",
        max_depth: Optional[int] = 1000,
        nonlin_coeff: Tuple[float, float] = (0, 0),
        signal_dict: Optional[SignalDict] = None,
        tstamp: Optional[float] = None,
        verbose: bool = False,
    ):
        """Initiate a SignalProcessor processor

        Args:
            mass_list (list of str): The names of the masses
            peak_type (str): The type of peak. Defaults to "gauss".
            max_depth (int or None): The max length of vectors in self.signal_dict
            nonlin_coeff (tuple of floats): The "nonlinear parameters" (P1, P2),
                which are the linear and quadratic dependence of sensitivity on p_vac
                in [10e-6 mbar]. Defaults to 0, i.e. linear. For parameters
                determined by Soren using closed-top chip measurements on LGACore in
                August 2020,  which approx. agree with those in Bela's linearity
                tests in March 2020, use the load constructor with 'Sorens final
                processor.json'
            signal_dict (SignalDict): The signal dictionary in which to record
                calculated signals
            tstamp (float): The unix timestamp which indicates the start of the measurement,
                used in the instantiation of signal_dict of not provided
            verbose (bool): Whether to print debugging info to terminal

        """
        self.mass_list = mass_list or []
        self.mass_peaks: Dict[str, Optional[Peak]] = {
            mass: None for mass in self.mass_list
        }
        self.signal_dict = signal_dict or SignalDict(
            max_depth=max_depth, tstamp=tstamp, verbose=verbose
        )
        self.peak_type = peak_type
        self.nonlin_coeff = nonlin_coeff
        self.medium = Medium()
        self.verbose = verbose

    @classmethod
    def load(
        cls,
        file_name: str = "Sorens final processor",
        proc_dir: Path = None,
        **kwargs: Any,
    ) -> "SignalProcessor":
        """Load the signal processor. The default has decent non-linearity correction."""
        proc_dir = proc_dir or CONFIG.processor_directory
        path_to_file = (proc_dir / file_name).with_suffix(".json")
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    @property
    def signals(self) -> SignalDict:
        """The SignalDict"""
        return self.signal_dict

    @property
    def p_vac(self) -> float:
        """The mediums p_vac"""
        return self.medium.p_vac

    @property
    def tstamp(self) -> float:
        """The unix start time of the measurement being processed"""
        return self.signal_dict.tstamp

    @property
    def PeakClass(self) -> Type[Peak]:
        """The type of Peak to use when calculating signals or other quantities"""
        if self.peak_type:
            return PEAK_CLASSES[self.peak_type]
        return Peak

    def calc_nonlinear_factor(self, p_vac: float = None) -> float:
        """Return the signal loss due to non-linearity given the vacuum pressure"""
        p_vac = p_vac if p_vac is not None else self.p_vac
        p_hat = p_vac / 1e-4  # pressure in 1e-4 Pa = 1e-6 mbar
        P1, P2 = self.nonlin_coeff
        nonlinear_factor = 1 + P1 * p_hat + P2 * p_hat**2
        if self.verbose:
            print(
                f"\tNonlinear factor = {nonlinear_factor} based on p_vac = {p_vac} [Pa]"
            )
        return nonlinear_factor

    def represent_nonlinear_correction(self) -> str:  # pragma: no cover
        """String representation of the non-linearity coefficients"""
        P1, P2 = self.nonlin_coeff
        return f"factor = 1 + {P1}(p_vac/[1e-6 mbar]) + {P2}(p_vac/[1e-6 mbar])^2"

    def __repr__(self) -> str:
        """The repr for this class"""
        self_as_str = f"{self.__class__}({self.represent_nonlinear_correction()})"
        return self_as_str

    def correct_y(
        self, x: FLOAT_ARRAY, y: FLOAT_ARRAY, p_vac: Optional[float] = None
    ) -> FLOAT_ARRAY:
        """Return the signal y corrected for background spectrum and non-linearity"""
        if self.verbose:
            print(f"SignalProcessor correcting signal in range m/z={[x[0], x[-1]]}.")
        nonlinear_factor = self.calc_nonlinear_factor(p_vac=p_vac)
        bg_interp = np.zeros(
            x.shape
        )  # Note: Implement here when background support is in
        y_corrected = (y - bg_interp) / nonlinear_factor
        return y_corrected

    def make_peak(
        self,
        x: FLOAT_ARRAY,
        y: FLOAT_ARRAY,
        mass: str = None,
        t: Optional[float] = None,
        Mspan: Optional[MSPAN] = None,
        fit_width: float = 1.0,
        p_vac: Optional[float] = None,
    ) -> Peak:
        """Return a background-subtracted peak for a given mass based on input data.

        The range of m/z values in the peak is either given by Mspan or determined
        by fit_width, centered at mass

        Side effects:
            Sets self.mass_peaks[mass]

        Args:
            x (np.array): m/z values for the raw data
            y (np.array): Signal values for the raw data in [A]
            mass (str): Mass for which to get the Peak, e.g. "M44"
            t (float): The measurement time with respect to self.tstamp of the peak.
                Defaults to now.
            Mspan (iterable): Optional. If given, [Mspan[0], Mspan[-1]] is taken to be
                the desired m/z range of the peak
            fit_width (float): The width of the peak in m/z units (if Mspan not given).
            p_vac (float): The vacuum pressure in [Pa]
        """
        if not Mspan:
            M = mass_to_M(mass)
            Mspan = (M - fit_width / 2, M + fit_width / 2)

        mask = np.logical_and(Mspan[0] <= x, x <= Mspan[-1])
        x = x[mask]
        y = y[mask]

        y = self.correct_y(x, y, p_vac=p_vac)
        t = t or time.time() - self.tstamp
        try:
            if self.verbose:
                print(f"fitting with {self.PeakClass}")
            peak = self.PeakClass(x=x, y=y, t=t)
        except PeakFitError as e:
            if self.verbose:
                print(
                    f"Failed to fit {self.PeakClass} at mass {mass} or mass span {Mspan} due to "
                    f"PeakFitError({e}). \n\tUsing simple Peak and recording error. "
                    f"Signal at this mass will be assumed to be zero."
                )
            peak = Peak(x=x, y=y, t=t, error=True)

        self.mass_peaks[mass] = peak
        return peak

    def calc_signal(
        self,
        x: FLOAT_ARRAY,
        y: FLOAT_ARRAY,
        mass: str,
        Mspan: Optional[MSPAN] = None,
        fit_width: float = 1.0,
        **kwargs: Any,
    ) -> float:
        """Calculate the advanced MID signal at mass as float in [A] based on raw data

        Uses ``self.make_peak(mass, x=x, y=y, Mspan=Mspan, width=width)`` if data given
        Otherwise uses the existing ``self.peaks[mass]``

        Side effects:

        * Sets ``self.peaks[mass]`` if data given
        * Calls ``self.peaks[mass].calc_signal()``, potentially fitting the peak.
        * Sets ``self.signal_dict[mass]``, thus appending to SignalDict vectors

        Args:
              x (np.array): The m/z values of the raw data
              y (np.array): The raw signal data in [A]
              mass (str): The mass string, e.g. "M44"
              Mspan (iterable): Optional. Custom mass range to be used for the peak
              fit_width (float): The width to use for make_peak()
              kwargs (dict): Additional key-word arguments are passed to
                  `Peak.calc_signal`

        Returns:
            float: The signal in [A]
        """
        if x is not None and y is not None:
            self.make_peak(x, y, mass=mass, Mspan=Mspan, fit_width=fit_width)
        S_M = self.mass_peaks[mass].calc_signal(**kwargs)
        self.signal_dict[mass] = S_M
        return S_M

    def calc_signals(
        self,
        x: FLOAT_ARRAY,
        y: FLOAT_ARRAY,
        mass_list: MASS_LIST,
        fit_width: float = 1.0,
    ) -> SignalDict:
        """Calculate and store a signal for each mass in mass_list based on raw data.

        Args:
            x (np.array): The m/z values of the raw data
            y (np.array): The raw signal data in [A]
            mass_list (list of str): The masses at which to calculate signal
            fit_width (float): The width to use for make_peak

        Returns:
             SignalDict: The signal dictionary containing the newly calculated signals
        """
        for mass in mass_list:
            self.calc_signal(x=x, y=y, mass=mass, fit_width=fit_width)
        return self.signal_dict

    def get_average_of_last(
        self, N: int, mass_list: MASS_LIST = None
    ) -> MASS_TO_SIGNAL:
        """Return ``{mass: S_M_avg}`` where s_M_avg is the average of mass's last N scans"""
        mass_list = mass_list or self.mass_list
        return self.signal_dict.get_average_of_last(N, mass_list)

    def __getattr__(self, attr: str) -> Dict[str, Any]:
        """Return the requested attr of the stored peak for each mass in mass_list"""
        try:
            if not self.mass_peaks:
                raise AttributeError
            return {mass: getattr(peak, attr) for mass, peak in self.mass_peaks.items()}
        except AttributeError:
            raise AttributeError(
                f"{self.__class__.__name__} either does not have any stored peaks or not all "
                f"of the peaks have the attribute '{attr}'"
            )
