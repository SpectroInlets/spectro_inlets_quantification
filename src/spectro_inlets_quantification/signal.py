# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Converting raw MS data to one S_vec per mass, corrected for bg and nonlinearity

The signal module is rather high up in the heiarchy so that it can be aware of
system conditiouns (Medium) and even the calibration (which I think is the right place
to store nonlinear coefficients.

The SignalProcesser, with help from the other classes of this and the peak module does:
    On initiation:
        * load one or more background spectra.
        * load one or more linearity correction functions
    During operation:
        * Process channel data to produce signal values (S_M) to feed to
            quant.calibration. This includes:
            - subtracting background from the channel data
            - fitting the resulting peaks
        * maybe adjusts backgrounds on the fly to match the peak-free areas of
            the measured full spectra?

The SignalDict and PeakSeries classes of this module are the (advanced) MID and spectrum
    data classes of quant. They have useful data visualization methods.

PeakSeries is the main analysis tool of recorded datasets (a parser should load raw data
    to a PeakSeries as the entry point to quant). It owns a SignalProcessor and uses it
    to generate the data in its SignalDict.

Likewise, a live data analyzer can use SignalProcessor to fill a SignalDict continuously
    during a measurement.
"""

from pathlib import Path  # noqa
import json
import numpy as np
from pandas import DataFrame
import time
from .constants import STANDARD_COLORS
from .tools import mass_to_M, mass_to_pure_mass, make_axis
from .exceptions import PeakFitError
from .medium import Medium
from .peak import PEAK_CLASSES, Peak
from .config import Config


CONFIG = Config()


class SignalDict:
    """Signals by (advanced) MID with memory. Index with mass for latest signal in [A].

    SignalDict is a wrapper around a set of mass-indexed signal vectors (S_vecs) and the
    corresponding time vectors (t_vecs) specifying when each signal was recorded. It
    is designed to be dynamic, making it easy to add MID signals as they are measured
    or calculated.

    Whenever a SignalDict is given a value assignment, it records the signal by
    appending it to S_vecs[mass] and appending the time to t_vecs[mass], removing
    the earliest one if needed due to the memory limitation set by max_depth.
    """

    def __init__(
        self,
        tstamp=None,
        max_depth=None,
        signals=None,
        verbose=False,
    ):
        """Initiate a SignalDict, optionally with initial signals

        Args:
            tstamp (float): The unix time defining the start of the measurement
            max_depth (int or None): The maximum length of the vectors, to save memory
            signals (dict): {mass: S_M_0} initial signals as float in [A] for each mass
            verbose (bool): Whether or not to print stuff to the terminal
        """
        self.tstamp = tstamp or time.time()
        self.S_vecs = {}
        self.t_vecs = {}
        self.max_depth = max_depth
        if signals:
            self.set_signals(signals, t=0)
        self.verbose = verbose

    @property
    def signals(self):
        """This is the dictionary with the latest signal for each mass"""
        S_dict = {mass: S_vec[-1] for mass, S_vec in self.S_vecs.items()}
        return S_dict

    def set_signal(self, mass, S_M, t=None):
        """Add the signal S_M at mass with measurement t (by default, now)"""
        t = t or (time.time() - self.tstamp)
        so_far = self.S_vecs.get(mass, np.array([]))
        t_so_far = self.t_vecs.get(mass, np.array([]))
        S_vec = np.append(so_far, S_M)
        t_vec = np.append(t_so_far, t)
        if self.max_depth and len(S_vec) > self.max_depth:
            S_vec = S_vec[-self.max_depth :]
            t_vec = t_vec[-self.max_depth :]
        self.S_vecs[mass] = S_vec
        self.t_vecs[mass] = t_vec

    def set_signals(self, signals, t=None):
        """Set multiple signals at once with the same measurement time"""
        for mass, S_M in signals.items():
            self.set_signal(mass, S_M, t=t)

    def get_signal(self, mass, tspan=None):
        """Return t, S_M, the time and signal vectors for mass during tspan"""
        t, S_M = self.t_vecs[mass], self.S_vecs[mass]
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t, S_M = t[mask], S_M[mask]
        return t, S_M

    def get_signals(self, mass_list=None, tspan=None):
        """Return the tspan-cut time and signal vecs for each mass in mass_list"""
        mass_list = mass_list or self.mass_list
        return {mass: self.get_signal(mass, tspan) for mass in mass_list}

    @property
    def mass_list(self):
        """The mass_list of a SignalDict lists the masses for which it has signals"""
        return list(self.S_vecs.keys())

    def signals_at(self, n):
        """Return the {mass: S_M} at the n'th iteration (negative is wrt most recent)"""
        return {mass: S_vec[n] for mass, S_vec in self.S_vecs.items()}

    def get_average_of_last(self, N=50, mass_list=None, t=None):
        """Return signal averaged over last N iterations or t seconds for mass_list"""
        mass_list = mass_list or self.mass_list
        if t:
            now = time.time() - self.tstamp
            tspan = [now - t, now]
            S_avg = {}
            for mass in mass_list:
                t_M, S_M = self.get_signal(mass, tspan)
                S_avg[mass] = np.mean(S_M)
        else:
            S_avg = {mass: np.mean(self.S_vecs[mass][-N:]) for mass in mass_list}
        return S_avg

    def __iter__(self):
        """Iterating a SignalDict iterates over the masses for which it has signals"""
        for mass in self.mass_list:
            yield mass

    def __contains__(self, mass):
        """A mass is in a SignalDict if it contains signal at that mass"""
        return mass in self.mass_list

    def items(self):
        """An iterator over (mass, most recent advanced MID signal at mass)"""
        for mass in self:
            yield mass, self[mass]

    def __repr__(self):
        return f"SignalDict({self.signals})"

    def __getitem__(self, key):
        """Index with a mass to get its most recent signal.

        Index with a negative integer to get the signals at a past iteration.
        Use cinf-style indexing for vectors: "{mass}-y" for signal and "{mass}-x" for
        the corresponding time vector.
        """
        if key in self:
            return self.signals[key]
        elif isinstance(key, int):
            return self.signals_at(key)
        elif key.endswith("-x"):
            return self.t_vecs[key[:-2]]
        elif key.endswith("-y"):
            return self.S_vecs[key[:-2]]
        raise KeyError

    def __setitem__(self, mass, value):
        """Set the signal with now as its measurement time"""
        self.set_signal(mass, value)

    def clear(self, mass):
        """Remove all knowledge of a mass from the SignalDict"""
        if mass in self:
            del self.S_vecs[mass]
            del self.t_vecs[mass]

    def clear_all(self):
        """Remove all knowledge at all masses from the SignalDict"""
        for mass in self:
            self.clear(mass)

    def plot(self, mass_list=None, tspan=None, ax="new"):
        """Plot the signal for each mass in mass_list"""
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("time / [s]")
            ax.set_ylabel("signal / [A]")
        mass_list = mass_list or self.mass_list
        for mass in mass_list:
            color = STANDARD_COLORS.get(mass_to_pure_mass(mass), None)
            t, S = self.get_signal(mass, tspan)
            ax.plot(t, S, color=color, label=mass)
        return ax


class SignalProcessor:
    """Class for calculating one signal (S_M) from the MS data for each peak (y vs x)

    The x-and-y stuff is actually handled by the Peak class and its family via make_peak
    Indexing a SignalProcessor with the mass M gives S_M.
    """

    def __init__(
        self,
        *,
        mass_list=None,
        peak_type="gauss",
        max_depth=1000,
        background=0,
        nonlin_coeff=(0, 0),
        medium=None,
        signal_dict=None,
        tstamp=None,
        verbose=False,
    ):
        """Initiate a SignalProcessor processor

        Args:
            background_file (str or Path): The file containing the x, y background data
                TODO: should have the possibility of different backgrounds for different
                    masses / pseudomasses
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
            medium (float): The Medium. medium.p_vac is vacuum chamber pressure in [Pa]
            signal_dict (SignalDict): the signal dictionary in which to record
                calculated signals
            verbose (bool): Whether to print debugging info to terminal.

        TODO: flesh out this class. Will require measurements and discussions.
        """
        self.background = background
        self.mass_list = mass_list or []
        self.mass_peaks = {mass: None for mass in self.mass_list}
        self.signal_dict = signal_dict or SignalDict(
            max_depth=max_depth, tstamp=tstamp, verbose=verbose
        )
        self.peak_type = peak_type
        self.nonlin_coeff = nonlin_coeff
        self.medium = medium or Medium()
        self.verbose = verbose

    @classmethod
    def load(cls, file_name="Sorens final processor", proc_dir=None, **kwargs):
        """Load the signal processor. The default has decent nonlinearity correction."""
        proc_dir = proc_dir or CONFIG.processor_directory
        path_to_file = (proc_dir / file_name).with_suffix(".json")
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    @property
    def signals(self):
        return self.signal_dict

    @property
    def p_vac(self):
        """SignalProcessor needs the Medium p_vac"""
        return self.medium.p_vac

    @property
    def tstamp(self):
        """The unix start time of the measurement being processed"""
        return self.signal_dict.tstamp

    @property
    def PeakClass(self):
        """The type of Peak to use when calculating signals or other quantities"""
        if self.peak_type:
            return PEAK_CLASSES[self.peak_type]
        return Peak

    def calc_nonlinear_factor(self, p_vac=None):
        """Return the signal loss due to nonlinearity given the vacuum pressure"""
        p_vac = p_vac if p_vac is not None else self.p_vac
        p_hat = p_vac / 1e-4  # pressure in 1e-4 Pa = 1e-6 mbar
        P1, P2 = self.nonlin_coeff
        nonlinear_factor = 1 + P1 * p_hat + P2 * p_hat**2
        if self.verbose:
            print(
                f"\tNonlinear factor = {nonlinear_factor} based on p_vac = {p_vac} [Pa]"
            )
        return nonlinear_factor

    def represent_nonlinear_correction(self):
        P1, P2 = self.nonlin_coeff
        return f"factor = 1 + {P1}(p_vac/[1e-6 mbar]) + {P2}(p_vac/[1e-6 mbar])^2"

    def __repr__(self):
        self_as_str = f"{self.__class__}({self.represent_nonlinear_correction()}, " + (
            "with background)" if self.has_background() else "without background)"
        )
        return self_as_str

    def has_background(self):
        if self.background is None:
            return False
        if type(self.background) in (int, float) and self.background == 0:
            return False
        return True

    def correct_y(self, x, y, p_vac=None):
        """Return the signal y corrected for background spectrum and nonlinearity"""
        if self.verbose:
            print(f"SignalProcessor correcting signal in range m/z={[x[0], x[-1]]}.")
        nonlinear_factor = self.calc_nonlinear_factor(p_vac=p_vac)
        bg_interp = np.zeros(x.shape)  # TODO: Actually get an interpolated background!
        y_corrected = (y - bg_interp) / nonlinear_factor
        return y_corrected

    def make_peak(
        self, x=None, y=None, t=None, mass=None, Mspan=None, fit_width=1, p_vac=None
    ):
        """Return a background-subtracted peak for a given mass based on input data.

            The range of m/z values in the peak is either given by mspan or determined
        by width, centered at mass

        Args:
            mass (str): mass for which to get the Peak, e.g. "M44"
            x (np.array): m/z values for the raw data
            y (np.array): signal values for the raw data in [A]
            t (float): The measurement time with respect to self.tstamp of the peak.
                Defaults to now.
            Mspan (iterable): optional. If given, [Mspan[0], Mspan[-1]] is taken to be
                the desired m/z range of the peak
            fit_width (float): the width of the peak in m/z units (if Mspan not given).
            p_vac (float): the vacuum pressure in [Pa]
        """
        if mass and not Mspan:
            M = mass_to_M(mass)
            Mspan = [M - fit_width / 2, M + fit_width / 2]
        if Mspan:
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
                    f"failed to fit {self.PeakClass} at {mass} due to "
                    f"PeakFitError({e}). \n\tUsing simple Peak and recording error. "
                    f"Signal at this mass will be assumed to be zero."
                )
                # raise  # debugging
            peak = Peak(x=x, y=y, t=t, error=True)
        self.mass_peaks[mass] = peak
        return peak

    def calc_signal(self, mass, x=None, y=None, Mspan=None, fit_width=1, **kwargs):
        """Calculate the advanced MID signal at mass as float in [A] based on raw data

        Uses self.make_peak(mass, x=x, y=y, Mspan=Mspan, width=width) if data given
        Otherwise uses the existing self.peaks[mass]

        Side effects:
            sets self.peaks[mass] if data given
            calls self.peaks[mass].calc_signal(), potentially fitting the peak.
            sets self.signal_dict[mass], thus appending to SignalDict vectors

        Args:
              mass (str): the mass string, e.g. "M44"
              x (np.array): the m/z values of the raw data
              y (np.array): the raw signal data in [A]
              Mspan (iterable): optional. Custom mass range to be used for the peak
              fit_width (float): the width to use for make_peak()
              kwargs (dict): additional key-word arguments are passed to

        Returns:
            float: The signal in [A]
        """
        if x is not None and y is not None:
            self.make_peak(mass=mass, x=x, y=y, Mspan=Mspan, fit_width=fit_width)
        S_M = self.mass_peaks[mass].calc_signal(**kwargs)
        self.signal_dict[mass] = S_M
        return S_M

    def calc_signals(self, mass_list, x=None, y=None, fit_width=1):
        """Calculate and store a signal for each mass in mass_list based on raw data.

        Args:
            mass_list (list of str): The masses at which to calculate signal
            x (np.array): the m/z values of the raw data
            y (np.array): the raw signal data in [A]
            fit_width (float): The width to use for make_peak

        Returns:
             SignalDict: The signal dictionary containing the newly calculated signals
        """
        for mass in mass_list:
            self.calc_signal(mass=mass, x=x, y=y, fit_width=fit_width)
        return self.signal_dict

    def adjust_background(self, full_scan_data):
        """TODO: discuss and implement this function"""
        pass

    def get_average_of_last(self, N, mass_list=None):
        """Return {mass: S_M_avg} where s_M_avg is the average of mass's last N scans"""
        mass_list = mass_list or self.mass_list
        return self.signal_dict.get_average_of_last(N, mass_list)

    def __getattr__(self, attr):
        """Return the requested attr of the stored peak for each mass in mass_list"""
        try:
            if not self.mass_peaks:
                raise AttributeError
            return {mass: getattr(peak, attr) for mass, peak in self.mass_peaks.items()}
        except AttributeError:
            raise AttributeError(
                f"SignalProcessor does not have attribute {attr} "
                f"and nor does it have mass_peaks which all have {attr}"
            )


# --------- the PeakSeries class -------------- #


class PeakSeries:
    """A class to track the evolution of a peak or a full spectrum over time."""

    def __init__(
        self,
        peak_list=None,
        df=None,
        x=None,
        mass=None,
        spectra=None,
        scalar_df=None,
        peak_type=None,
        signal_processor=None,
        tstamp=None,
        verbose=False,
    ):
        """Initiate a peak_series from either a peak_list or a DataFrame

            It can be initiated from a series of spectra (peak_list), or from a
        dataframe (df) with a 2-d array containing the spectrum data, together with the
        corresponding m/z values (x).

        The mass argument is important. If provided, then the PeakSeries is take to
        track a single mass, opening up the possibility to fit each peak in the series
        without further sniping in m/z. Otherwise,

        Args:
            peak_list (list of Peak): list of the peaks in the series
            df (pandas.DataFrame): the intensity data from the spectra
                (if peak_list not given)
            spectra (np.array): the intensity data from the spectra (if peak_list and
                df not given)
            x (numpy.array): m/z-values corresponding to the peaks
            scalar_df (pandas.DataFrame): scalar data corresponding to the peaks
                This is, e.g. time "t", or associated data like "pressure".
                If not given, scalar_df will just contain a trivial "counter"
            peak_type (str): default peak type, only used if no signal_processor given
            signal_processor (SignalProcessor): The signal processor used to evt.
                subtract background and correct for nonlinearity
            tstamp (float): The unix timestamp, by default use the signal_processor's
        """
        self.signal_processor = signal_processor or SignalProcessor(
            tstamp=tstamp, verbose=verbose, peak_type=peak_type
        )
        self.signal_dict = self.signal_processor.signal_dict
        if tstamp:
            self.signal_dict.tstamp = tstamp
        self._mass = mass
        self.x = x
        if peak_list:
            spectra = None
            # generate the 2-D array with data for all the spectra stacked
            for peak in peak_list:
                y = peak.y
                # stack the spectra:
                spectra = np.append(spectra, y) if spectra else y
            df = DataFrame(spectra)
            if not self.x:
                self.x = peak_list[0].x

        self._peak_list = peak_list
        if df is None and spectra is not None:
            df = DataFrame(spectra)
        self.df = df
        if scalar_df is None:
            scalar_df = DataFrame({"counter": np.arange(len(self.df))})
        self.scalar_df = scalar_df
        self.S_vecs = {}  # to store mass-specific signals (for full-scan channels)
        self._S_vec = None  # to store full-channel signals (e.g. single-peak channels)
        self._making_peak = False  # prevent recursion between __getattr__() and peaks()
        self.verbose = verbose

    @property
    def spectra(self):
        return self.df.to_numpy()

    @property
    def tstamp(self):
        """The unix time of the measurement start"""
        return self.signal_dict.tstamp

    @property
    def peak_type(self):
        """The type of peak to use when calculating signals or other quantities"""
        return self.signal_processor.peak_type

    @property
    def PeakClass(self):
        """The Peak class to use when calculating signals or other quantities"""
        return self.signal_processor.PeakClass

    def get_peak(self, n, Mspan=None, mass=None, fit_width=1.0):
        """Return a corrected and fitted peak centered at mass for the n'th scan

        If given a Mspan or mass, it cuts the data in m/z before fitting a peak.
        Otherwise, it fits the whole m/z-range, only meaningful for single-mass scans.

        FIXME: This method does not check that it is a single-mass peak. The code
            # mass = mass or self.mass
            would do it, but that cause it to fit over fit_width centered on mass
            rather than the "natural" Mspan built into the data which is prefered for a
            true PeakSeries of single-mass scans
        """
        self._making_peak = True  # used to avoid an infinite recursion.
        y = self.df.to_numpy()[n]
        try:
            t = self.t[n]
        except AttributeError:
            t = None
        try:
            p_vac = self.p_vac[n]
        except AttributeError:
            p_vac = None
        #  p_vac = self.scalar_df.get("p_vac", {key: None})[key]  # <- alternative
        peak = self.signal_processor.make_peak(
            x=self.x, y=y, p_vac=p_vac, t=t, Mspan=Mspan, mass=mass, fit_width=fit_width
        )
        self._making_peak = False
        return peak

    def __getitem__(self, key):
        """Return a single peak (int index) or a scalar from all peaks (str index)"""
        if isinstance(key, int):
            return self.get_peak(key)
        try:
            return self.signal_dict[key]
        except KeyError:
            try:
                return self.scalar_df[key].to_numpy()
            except KeyError:
                raise KeyError(
                    f"The key '{key}' can be found in neither the signal_dict nor the "
                    f"scalar_df of the PeakSeries"
                )

    def __len__(self):
        """The number of scans in the PeakSeries"""
        return len(self.df)

    def peaks(self):
        """Iterate over the represented peaks. Meaningful for single-mass PeakSeries"""
        for n in range(len(self)):
            yield self.get_peak(n)

    def __getattr__(self, attr):
        """Look for attr in scalar_df, otherwise get attr from each peak for vector.

        Unnamed attributes are taken to be numerical metadata associated with each peak
        in the series. These are generally stored in scalar_df.
        If attr is not in self.scalar_df, try and get the attr from each peak and
        cache the resulting vector as self.scalar_df[attr] before returning it.
        """
        if attr in self.scalar_df:
            if attr == "counter":
                return self.scalar_df[attr].to_numpy(dtype=int)
            return self.scalar_df[attr].to_numpy()
        try:
            if self._making_peak:
                raise AttributeError  # this is to avoid an infinite recursion
            values = []
            for peak in self.peaks():
                value = getattr(peak, attr)
                if attr == "signal" and peak.error:
                    value = 0
                values += [value]
        except AttributeError:
            raise AttributeError(
                f"PeakSeries couldn't get {attr} from its scalar_df or peaks()"
            )
        else:
            values = np.array(values)
            self.scalar_df[attr] = values  # cache the result
        return values

    @property
    def mass(self):
        if not self._mass:
            raise TypeError(
                "PeakSeries method wrongly called which only works on a single-mass"
                "peak series, i.e. one in which mass was set on initiation"
            )
        return self._mass

    @property
    def S_vec(self):
        """The cache'd or calculated signal vector, useful for single-peak PeakSeries"""
        mass = self.mass
        if mass not in self.signal_dict:
            self.calc_signal()
        return self.signal_dict[mass + "-y"]

    def set_peak_type(self, peak_type):
        """Changing the peak type clears the signal dict."""
        self.signal_processor.peak_type = peak_type
        self.signal_dict.clear_all()

    def cut(self, tspan):
        """Return a PeakSeries cut according to the tspan"""
        t = self.t
        mask = np.logical_and(tspan[0] < t, t < tspan[-1])
        df_cut = self.df[mask]
        df_scalars_cut = self.scalar_df[mask]
        return PeakSeries(
            df=df_cut,
            x=self.x,
            scalar_df=df_scalars_cut,
            signal_processor=self.signal_processor,
        )

    def snip(self, mass=None, fit_width=1.0, Mspan=None):
        """Return a PeakSeries cross section according to an m/z range and peak_type"""
        if mass and not Mspan:
            M = mass_to_M(mass)
            Mspan = [M - fit_width / 2, M + fit_width / 2]
        if not Mspan:
            raise TypeError("PeakSeries.snip requires either mass or Mspan")
        mask = np.logical_and(Mspan[0] < self.x, self.x < Mspan[-1])
        snipped_x = self.x[mask]
        snipped_df = self.df.iloc[:, mask]
        return PeakSeries(
            df=snipped_df,
            x=snipped_x,
            mass=mass,
            scalar_df=self.scalar_df.copy(),
            signal_processor=self.signal_processor,
            verbose=self.verbose,
        )

    def calc_signal(self, mass=None, fit_width=1.0, peak_type=None):
        """Calculate, cache, and return the np.array for the signal at mass

        Args:
            mass (str): the mass at which to calculate the signal vector. Not needed
                for single-mass PeakSeries
            fit_width (float): If mass is given, the width of the m/z range centered on
                mass to include in the fit
            peak_type (str): The fitting type. Default: self.signal_processor.peak_type

        Side-effects:
            Cache's the calculated vector. If mass is given, it is cached as
                self.signals[mass]. If mass is not given (single-peak PeakSeries), it
                is cache'd as self._S_vec, accessible as property self.signal.
        """
        if mass:
            peak_series = self.snip(mass=mass, fit_width=fit_width)
        else:
            mass = self.mass  # only available for single-mass PeakSeries
            peak_series = self
        if peak_type:
            peak_series.signal_processor.peak_type = peak_type
        self.signal_dict.clear(mass)  # removes anything that was already there.
        for peak in peak_series.peaks():  # already processed by self.signal_processor!
            S_n = peak.calc_signal()  # will give zero if the peak could not be fit!
            self.signal_dict.set_signal(mass, S_n, t=peak.t)
        return self.signal_dict.get_signal(mass)

    def get_signal(self, mass=None, tspan=None, fit_width=1.0, peak_type=None):
        """Return np.arrays for t and signal at mass, cache'd else calculated"""
        if not mass:
            mass = self.mass
        if mass not in self.signal_dict:
            self.calc_signal(mass, fit_width=fit_width, peak_type=peak_type)
        return self.signal_dict.get_signal(mass, tspan=tspan)

    def get_signals(self, mass_list, tspan=None, fit_width=1.0, peak_type=None):
        """Return dict of {mass: (t, S_M)} where S_M is the signal at mass"""
        signals = {}
        for mass in mass_list:
            signals[mass] = self.get_signal(
                mass, tspan=tspan, fit_width=fit_width, peak_type=peak_type
            )
        return signals

    def plot_signals(
        self,
        mass_list,
        tspan=None,
        fit_width=1.0,
        vs="t",
        peak_type=None,
        ax=None,
        logscale=False,
        legend=False,
        **kwargs,
    ):
        """Plot the series of spectra in self.df as advanced MID

        Uses cache'd signals if available else calculates them by fitting peaks

        Args:
            mass_list (list of str): The masses to plot
            tspan (iterable): The time range to plot
            fit_width (float): The width to assume if fitting peaks
            vs (str): The attribute to put on the x-axis
            peak_type (str): The type of peak fitting to use if calculating peaks
                By default, it's self.signal_processer.peak_type
            ax (matplotlib axis): The axis to plot on. Makes a new axis by default
            logscale (bool): Whether to plot on a log scale
            legend (bool): Whether to include a legend
            kwargs (dict): Additional key-word arguments go to matplotlib.pyplot.plot()
        """
        if not ax:
            fig, ax = make_axis()
            ax.set_xlabel("time / [s]" if vs == "t" else vs)
            ax.set_ylabel("signal / [A]")
        t_substitute = getattr(self, vs)
        signals = self.get_signals(mass_list, tspan, fit_width, peak_type)
        for mass, (t, S_vec) in signals.items():
            t_like = t if vs == "t" else t_substitute
            color = STANDARD_COLORS.get(mass_to_pure_mass(mass), None)
            ax.plot(t_like, S_vec, color=color, label=mass, **kwargs)
        if logscale:
            ax.set_yscale("log")
        if legend:
            ax.legend()
        return ax

    def plot_spectrum(
        self, mode="average", logscale=False, color="k", ax="new", shift=None
    ):
        """plot the peak series as a single spectrum."""
        # if not mode == "average":
        #    raise NotImplementedError
        # print("jj")

        if ax == "new":
            fig, ax = make_axis()
        spectra = self.spectra

        if mode == "average":
            y_avg = spectra.mean(axis=0)

            ax.plot(self.x, y_avg, color=color)
        elif mode == "stacked":
            for n, y in enumerate(spectra):
                if logscale:
                    ax.plot(self.x, y * shift**n)
                else:
                    ax.plot(self.x, y - shift * n)
        else:
            raise NotImplementedError(f"Mode '{mode}' not implemented")

        if logscale:
            ax.set_yscale("log")
        return ax

    def heat_plot(
        self, orientation="yx", vs=None, logscale=False, zrange=None, **kwargs
    ):
        """Plot the peak series as a heat map

        Args:
            orientation (str): "yx" puts m/z on x-axis. "xy" puts m/z on y-axis.
            vs (str): What to put on second axis of heat plot (first being m/z). Must
                be attribute-accessible - for example "t" or "counter". "t" by default.
            logscale: whether to take log of intensity before applying the color scale
            zrange: range outside of which MS data is removed. By default only negatives
                and np.nan's are removed.
            kwargs: other key-word arguments are passed to matplotlib.pyplot.imshow()
        """

        spectra = self.spectra
        x = self.x

        if not vs:
            vs = "t" if "t" in self.scalar_df else "counter"
        t_like = getattr(self, vs)
        t_label = "time / [s]" if vs == "t" else vs

        fig, ax = make_axis()

        if orientation == "xy":
            spectra = np.swapaxes(spectra, 0, 1)
            spectra = np.flip(spectra, axis=0)
            ax.set_ylabel("m/z")
            ax.set_xlabel(t_label)
        else:
            ax.set_ylabel(t_label)
            ax.set_xlabel("m/z")
        if "extent" not in kwargs:
            if orientation == "xy":
                kwargs["extent"] = [t_like[0], 2 * t_like[-1] - t_like[-2], x[0], x[-1]]
            elif orientation == "yx":
                kwargs["extent"] = [x[0], x[-1], 2 * t_like[-1] - t_like[-2], t_like[0]]

        if logscale:
            spectra = np.log(spectra)
        if zrange is None:
            good = np.logical_and(~np.isnan(spectra), ~np.isinf(spectra))
            # print('spectra = \n' + str(spectra)) # debugging
            low = np.min(spectra[good])
            high = np.max(spectra[good])
        else:
            low = zrange[0]
            high = zrange[1]

        spectra[spectra < low] = low
        spectra[spectra > high] = high
        spectra[np.isnan(spectra)] = low
        spectra[np.isinf(spectra)] = low

        if "aspect" not in kwargs:
            kwargs.update(aspect="auto")
        elif kwargs["aspect"] == "square":
            if orientation == "xy":
                kwargs["aspect"] = (t_like[-1] - t_like[0]) / (x[-1] - x[0])
            elif orientation == "yx":
                kwargs["aspect"] = (x[-1] - x[0]) / (t_like[-1] - t_like[0])
        if "cmap" not in kwargs:
            kwargs.update(cmap="inferno")

        ax.imshow(spectra, **kwargs)
        return ax
