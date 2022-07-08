# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""The SensitivityFactor, SensitivityMatrix, and QuantificationMatrix classes

A few abbreviations used throughout this module:
    mdict - molecule dictionary, the one and only instance of MoleculeDict
    sf - sensitivity factor, an instance of SensitivityFactor
    sf_list - a list of SensitivityFactor instances
    sensitivity_list - sensitivity list, an instance of SensitivityList
    sm - sensitivity matrix, an instance of SensitivityMatrix
    fit - sensitivity fit, an instance of SensitivityFit
    cal - calibration point, an instance of CalPoint
    cal_list - a list of CalPoint instances

F is absolute sensitivity (signal / molecular flux) in [A / (mol/s)] = [C/mol]
f is predicted sensitivity relative to the sensitivity for O2 at m/z=32

See MS_Theory_v1p0.pdf for details.

FIXME Insert reference to quantitative MS docs

"""

from math import isclose

import attr
import numpy as np
from scipy.optimize import minimize
from .constants import (
    REFERENCE_MASS,
    REFERENCE_MOLECULE,
    STANDARD_TRANSMISSION_EXPONENT,
    STANDARD_IONIZATION_ENERGY,
    CAL_TYPE_SPECS,
    STANDARD_COLORS,
    STANDARD_MOL_COLORS,
)
from .tools import (
    make_axis,
    mass_to_M,
    mass_to_setting,
    mass_to_pure_mass,
    dict_equal_with_close_floats,
)
from .exceptions import SensitivityMatrixError
from .molecule import MoleculeDict


@attr.s
class SensitivityFactor:
    """Data class for storing an f and F

    Attributes:
        mol (str): name of molecule
        mass (str): name of mass
        F (float): F^{mol}_{mass}, the absolute sensitivity factor in [C/mol]
        f (float): The predicted relative sensitivity factor, used by SensitivityFit
        F_type (str): Description of the origin of the sensitivity factor, e.g.
            "internal", "external", "predicted", etc.
    """

    mol = attr.ib()
    mass = attr.ib()
    F = attr.ib()
    f = attr.ib(default=None)
    F_type = attr.ib(default=None)

    @property
    def pure_mass(self):
        return mass_to_pure_mass(self.mass)

    @property
    def setting(self):
        return mass_to_setting(self.mass)

    @property
    def M(self):
        return mass_to_M(self.mass)

    def union(self, other):
        """Return self's SensitivityUnion with a SensitivityUnion or SensitivityFactor

        self and other need to have the same mol and mass.
        """
        if not (self.mol == other.mol and self.mass == other.mass):
            raise TypeError(f"can't combine {self} and {other}")
        if isinstance(other, SensitivityUnion):
            return other.union(self)
        return SensitivityUnion(
            mol=self.mol, mass=self.mass, sf_list=[self, other], f=self.f
        )

    def as_dict(self):
        """Return the dictionary representation of self"""
        self_as_dict = dict(
            mol=self.mol, mass=self.mass, F=self.F, f=self.f, F_type=self.F_type
        )
        return self_as_dict

    def copy(self):
        """Return a new SenitivityFactor (or inheriting) object cloning self"""
        return SensitivityFactor(
            mass=self.mass, mol=self.mol, F=self.F, f=self.f, F_type=self.F_type
        )


class SensitivityUnion(SensitivityFactor):
    """A class combining multiple SensitivityFactors of the same mol and mass

    The idea behind this is to make it easy to keep track of the accuracy of each
    sensitivity factor, i.e. the variation of the same sensitivity factor measured in
    different ways. The original SensitivityFactors are listed here in sf_list.

    A sensitivity union can be used just like a normal sensitivity factor.
    """

    def __init__(self, mol, mass, sf_list, f):
        """Initiate the sensitivity union

        self.F becomes the average F of the sensitivity factors in sensitivity_list
        self.F_error (the accuracy) becomes the standard deviation thereof.
        self.F_type is "union"

        args:
            mol (str): name of molecule
            mass (str): name of mass
            sensitivity_list (list of SensitivityFactors): the individual SensitivityFactors
            f (float): The predicted relative sensitivity factor, used by SensitivityFit
        """
        F_vec = np.array([sf.F for sf in sf_list])
        F = np.mean(F_vec)
        super().__init__(mol=mol, mass=mass, F=F, f=f, F_type="union")
        self.sf_list = sf_list
        self.F_error = np.std(F_vec)

    def as_dict(self):
        """Return the dictionary representation, with the individual sf's"""
        self_as_dict = super().as_dict()
        sf_dicts = [sf.as_dict() for sf in self.sf_list]
        self_as_dict.update(sf_list=sf_dicts, F_error=self.F_error)
        return self_as_dict

    def union(self, other):
        """Return self's SensitivityUnion with a SensitivityUnion or SensitivityFactor

        self and other need to have the same mol and mass.
        """
        if not self.mol == other.mol and self.mass == other.mass:
            raise TypeError(f"can't combine {self} and {other}")
        if isinstance(other, SensitivityUnion):
            sf_list = self.sf_list + other.sf_list
        elif isinstance(other, SensitivityFactor):
            sf_list = self.sf_list + [other]
        else:
            raise TypeError(f"can't combine {self} and {other}")
        return SensitivityUnion(mol=self.mol, mass=self.mass, sf_list=sf_list, f=self.f)

    def copy(self):
        """Return a cloned SensitivityUnion of cloned SensitivityFactors"""
        new_sf_list = [sf.copy() for sf in self.sf_list]
        return SensitivityUnion(
            mol=self.mol, mass=self.mass, sf_list=new_sf_list, f=self.f
        )

    @property
    def accuracy(self):
        """The relative standard deviation of the contained sensitivity factors"""
        return self.F_error / self.F


class SensitivityList:
    """A wrapper around a list of SensitivityFactors of various mol & mass"""

    def __init__(self, sf_list):
        """Initiate a SensitivityList from a list of SensitivityFactors"""
        self.sf_list = sf_list

    def __getitem__(self, key):
        return self.sf_list[key]

    def __len__(self):
        return len(self.sf_list)

    def __add__(self, other):
        """Return a SensitivityList with the sensitivity factors of self and other"""
        if isinstance(other, SensitivityList):
            sf_list = self.sf_list + other.sf_list
        else:
            raise TypeError(f"can't add {other} to {self}")
        return SensitivityList(sf_list)

    def __iadd__(self, other):
        """Add the SensitivityFactors from a SensitivityList or list"""
        if isinstance(other, SensitivityList):
            self.sf_list = self.sf_list + other.sf_list
        elif isinstance(other, (list, tuple)):
            for sf in other:
                self.append(sf)
        else:
            raise TypeError(f"can't add {other} to {self}")
        return self

    def append(self, sf):
        """Add a SensitivityFactors to self's list"""
        if not isinstance(sf, SensitivityFactor):
            raise TypeError(
                f"Can only append SensitivityFactor instances to "
                f"{self.__class__}. You tried to append {sf}."
            )
        self.sf_list.append(sf)
        return self

    def __iter__(self):
        yield from self.sf_list

    def __repr__(self):
        """A somewhat verbose string representation showing all the contents"""
        self_as_str = f"{self.__class__}(["
        for sf in self:
            self_as_str += f"\n\t{sf},"
        self_as_str += "\n])"
        return self_as_str

    def to_sf_dict(self):
        """Return a dictionary rearranging the sensitivity factors by [mol][mass]

        This also joins duplicate sensitivity factors by their union() method.
        """
        sf_dict = {}
        for sf in self.sf_list:
            mol = sf.mol
            mass = sf.mass
            if mol not in sf_dict:
                sf_dict[mol] = {}
            if mass not in sf_dict[mol]:
                sf_dict[mol][mass] = sf
            else:
                sf_dict[mol][mass] = sf_dict[mol][mass].union(sf)
        return sf_dict

    def to_sensitivity_matrix(self, mol_list, mass_list, fit=None, metadata=None):
        """Return a SensitivityMatrix instance based on the contained SensitivityFactors

        Args:
            mol_list (list of str): the molecules of the SensitivityMatrix
            mass_list (list of str): the masses of the SensitivityMatrix
            fit (SensitivityFit): The fit to use when predicting missing F's
            metadata (dict): the metadata of the SensitivityMatrix
        """
        full_sf_dict = self.to_sf_dict()
        sf_dict = {
            mol: {mass: sf for mass, sf in full_sf_dict_i.items() if mass in mass_list}
            for mol, full_sf_dict_i in full_sf_dict.items()
            if mol in mol_list
        }
        return SensitivityMatrix(
            sf_dict=sf_dict,
            mol_list=mol_list,
            mass_list=mass_list,
            fit=fit,
            metadata=metadata,
        )

    def filter(self, **kwargs):
        """Return a SensitivityList with a subset of the contained SensitivityFactors.

        Best described with an example:
        >>> sensitivity_list = SensitivityList([
        ...     SensitivityFactor(mol="H2", mass="M2", F=4, F_type="internal"),
        ...     SensitivityFactor(mol="CO2", mass="M28", F=2, F_type="semi"),
        ...     SensitivityFactor(mol="ethanol", mass="M43", F=3, F_type="predicted")
        ...])
        >>> sensitivity_list.filter(mol="CO2")
        SensitivityList([SensitivityFactor(mol="CO2", mass="M28", F=2, F_type="semi")])
        >>>

        A string argument starting with a "!" filters against rather than for the value:
        >>> sensitivity_list.filter(F_type="!predicted")
        SensitivityList([
            SensitivityFactor(mol="H2", mass="M2", F=4, F_type="internal"),
            SensitivityFactor(mol="CO2", mass="M28", F=2, F_type="semi"),
        ])
        >>>
        """
        new_sf_list = []
        for sf in self:
            good = True
            for key, values in kwargs.items():
                if isinstance(values, str) and values.startswith("!"):
                    good = good and (getattr(sf, key) != values[1:])
                else:
                    if not isinstance(values, (list, tuple)):
                        values = [values]
                    if not getattr(sf, key) in values:
                        good = False
            if good:
                new_sf_list.append(sf)
        return SensitivityList(new_sf_list)

    def as_dict(self):
        """Return a full dictionary representation of the sensitivity list"""
        sf_dicts = [sf.as_dict for sf in self.sf_list]
        self_as_dict = {"sensitivity_list": sf_dicts}
        return self_as_dict


class SensitivityMatrix:
    """Class for handling and using an array of sensitivity factors for quantification

    SensitivityMatrix.F_mat is the sensitivity factor matrix itself as a numpy array.
    SensitivityMatrix acts like F_mat in that indexing with integers and slices
    returns the corresponding number, row, or column, in F_mat (example below).

    SensitivityMatrix.Q_mat is the quantification matrix, the (pseudo)inverse to F_mat

    Indexing with strings interprets the first index as the molecule name (mol) and
    the second as the mass.
    SensitivityMatrix.molecule() returns a molecule "calibrated" with its own
    quantification vector that it can dot with a corresponding S_vec vector to return
    its flux (will be very useful for making calibrated EC-MS plots).

    Normally a sensitivity matrix will be generated not directly, but by
        SensitivityList.to_sensitivity_matrix

    Otherwise, these examples show it's use:
    >>> sm = SensitivityMatrix(
    ...     mol_list=["H2", "O2"], mass_list=["M2", "M32"],
    ...     sf_dict={
    ...         "H2": {"M2": SensitivityFactor("H2", "M2", 1)},
    ...         "O2": {"M32": SensitivityFactor("O2", "M32", 2)}
    ...     }
    ... )
    >>> sm[0]
    array([1., 0.])
    >>> sm[1]
    array([0., 2.])
    >>> sm["H2"]
    {'M2': SensitivityFactor(mol='H2', mass='M2', F=1, f=None, F_type=None),
    'M32': SensitivityFactor(mol='H2', mass='M32', F=0.0, f=0, F_type='predicted')}
    >>> sm["O2"]
    {'M32': SensitivityFactor(mol='O2', mass='M32', F=2, f=None, F_type=None),
     'M2': SensitivityFactor(mol='O2', mass='M2', F=0.0, f=0, F_type='predicted')}
    >>> sm.Q_mat
    array([
        [1. , 0. ],
        [0. , 0.5]
    ])
    >>>
    """

    def __init__(
        self,
        *,
        mol_list=None,
        mass_list=None,
        sf_dict=None,
        fit=None,
        fit_specs=None,
        metadata=None,
        verbose=False,
    ):
        """Initiate a SensitivityMatrix

        Args:
            mol_list (list of str): The molecules of the sensitivity matrix
                These corresponds to the columns of F_mat and rows of Q_mat
            mass_list (list of str): The masses of the sensitivity matrix
                These correspond to the rows of F_mat and the columns of Q_mat
            sf_dict (dict of dicts of SensitivityFactors): two-layer dictionary with
                sf_dict[mol][mass] = SensitivityFactor(mol=mol, mass=mass, ...)
            fit (SensitivityFit): the fit, used to predict any missing F's
            fit_spec (dict): The parameters to pass into SensitivityFit to create in the fit,
                instead of using the one provided as `fit`. This consists of the `as_dict`
                return value.
            metadata (dict): associated metadata, used by the Calibration class.
            verbose (bool): Whether to print stuff to terminal
        """
        # TODO Rename fit_spec to follow similar naming elsewhere on the form fit_as_dict
        self.mol_list = mol_list
        self.mass_list = mass_list
        self.sf_dict = sf_dict
        if fit_specs:
            self.fit = SensitivityFit(**fit_specs)
        else:
            self.fit = fit
        self.metadata = metadata
        self.mdict = MoleculeDict()
        self._Q_mat = None  # to be calculated later
        self.verbose = verbose

    def __eq__(self, other):
        """Returns whether a SensitivityMatrix is equal to another"""
        if self.fit != other.fit:
            return False
        if not dict_equal_with_close_floats(self.sf_dict, other.sf_dict):
            return False
        return True

    def __repr__(self):
        return (
            f"SensitivityMatrix(mol_list={self.mol_list}, mass_list={self.mass_list})"
        )

    def __getitem__(self, key):
        """Indexing returns from the sensitivity matrix or sensitivity factor stack"""
        if isinstance(key, (int, slice, tuple)):
            return self.F_mat[key]
        elif isinstance(key, str):
            return self.sf_dict[key]
        raise KeyError(
            f"Key to SensitivityMatrix must be int, slice, tuple, or str. Got {key}"
        )

    def as_dict(self):
        """Return the dictionary representation of self

        NOTE: Saving and loading is handled by Calibration, but that needs a dict

        NOTE: This is, as of yet, not JSON-able, but is picklable.
        """
        # TODO: For the full JSON-able implementation, have a look here:
        #  https://github.com/SpectroInlets/spitze/pull/84#discussion_r889043589

        # sf_as_dict_list = [sf.as_dict() for sf in self.to_sensitivity_list().sf_list]
        self_as_dict = dict(
            # sf_as_dict_list=sf_as_dict_list,
            sf_dict=self.sf_dict,
            # fit_specs becomes dict with {'sensitivity_list' -> [sf_dict, ...]}
            fit_specs=self.fit.as_dict(),
            metadata=self.metadata,
        )
        return self_as_dict

    def to_sensitivity_list(self):
        """Return a SensitivityList with all of the SensitivityFactors in self"""
        sl = SensitivityList([])
        for mol in self.sf_dict:
            for mass in self.sf_dict[mol]:
                sl.append(self.sf_dict[mol][mass])
        return sl

    @property
    def N_mol(self):
        """The number of molecules in the sensitivity matrix"""
        return len(self.mol_list)

    @property
    def N_mass(self):
        """The number of masses in the sensitivity matrix"""
        return len(self.mass_list)

    @property
    def F_mat(self):
        """The sensitivity matrix. Signal per flux in [C/mol]. Built each time."""
        return np.array(
            [
                [self.get_F(mol, mass) for mol in self.mol_list]
                for mass in self.mass_list
            ]
        )

    def prints_F_mat(self):
        """Format the active sensitivity factor matrix and return the string"""
        return str(np.round(self.F_mat, decimals=3))

    def print_F_mat(self):
        """Print the active sensitivity factor matrix in decent-looking format."""
        print(self.prints_F_mat())

    def make_fit(self, **kwargs):
        sf_fit = (
            self.to_sensitivity_list()
            .filter(**kwargs)
            .filter(F_type="!predicted")
            .filter(F_type="!ref_spectrum")
        )
        self._fit = SensitivityFit(sf_fit)
        self._fit.fit()

    @property
    def fit(self):
        """The sensitivity matrix's fit of F vs f for its sensitivity factors"""
        if not self._fit:
            self.make_fit()
        return self._fit

    @fit.setter
    def fit(self, new_fit):
        """Recalibration requires re-setting the fit"""
        self._fit = new_fit

    @property
    def alpha(self):
        """The overall sensitivity number"""
        return self.fit.alpha

    @property
    def beta(self):
        """The exponent of the mass-dependence of sensitivity factors"""
        return self.fit.beta

    def fit_F_vs_f(self, **kwargs):
        """fit the SensitivityFactor trend with self.fit"""
        if not self.fit:
            self.make_fit(**kwargs)
        self.fit.fit()

    def plot_F_vs_f(self, **kwargs):
        """fit the SensitivityFactor trend with self.fit"""
        if not self.fit:
            self.make_fit(**kwargs)
        return self.fit.plot_F_vs_f()

    def get_F(self, mol, mass):
        """Return the most trusted available F^{mol}_{mass}.

        Priority:
            (1) saved already in the matrix, i.e. from self.sf_dict[mol][mass]
            (2) a saved factor for another mass in the molecule, scaled according to
                the molecule's spectrum
            (3) a sensitivity factor predicted by self.fit

        Args:
            mol (str): molecule name
            mass (str): mass

        Returns float: sensitivity factor in [C/mol]
        """
        my_spectrum = None
        if mol in self.sf_dict:
            if mass in self.sf_dict[mol]:
                # beautiful! we already have the sensitivity factor. Return it.
                return self.sf_dict[mol][mass].F
            # okay, we don't have the exact sensitivity factor, but we have another
            #   for the molecule. Lets assume its reference spectrum is right and
            #   scale the (biggest) one we have accordingly.
            my_spectrum = {
                mass: sf.F
                for mass, sf in self.sf_dict[mol].items()
                if "predicted" not in sf.F_type
            }
        if my_spectrum:  # Then there is a measured value for mol at another mass
            ref_spectrum = self.mdict.get(mol).norm_spectrum
            # We're only interested in masses at the same settings, otherwise we
            #   don't trust that the spectrum translates. Here we filter:
            try:
                masses, Fs = zip(
                    *[
                        (mass_, F)
                        for mass_, F in my_spectrum.items()
                        if mass_to_setting(mass_) == mass_to_setting(mass)
                    ]
                )
            except ValueError:
                if self.verbose:
                    print(
                        f"Could not find appropriate reference values in my_spectrum = "
                        f"{my_spectrum} to use ref_spectrum = {ref_spectrum} "
                        f"to calculate F({mol}, {mass})."
                    )
            else:
                # go for the closest mass:
                I_closest = np.argmin(
                    [np.abs(mass_to_M(mass_) - mass_to_M(mass)) for mass_ in masses]
                )
                mass_closest = masses[I_closest]
                F_closest_mass = Fs[I_closest]
                if mass in ref_spectrum and mass_closest in ref_spectrum:
                    F = F_closest_mass * ref_spectrum[mass] / ref_spectrum[mass_closest]
                    sf = SensitivityFactor(
                        mass=mass, mol=mol, F=F, f=None, F_type="ref_spectrum"
                    )
                    self.sf_dict[mol][mass] = sf
                    if self.verbose:
                        print(f"got {sf} via reference spectrum!")
                    return F
        # If we get here it means that there's no saved sensitivity data on [mol], or
        #   there is but it can't be used with the reference spectrum for mol. So we
        #   have to calculate it with the fit.
        sf = self.fit.predict_sf(mol, mass)
        if mol not in self.sf_dict:
            self.sf_dict[mol] = {}
        self.sf_dict[mol][mass] = sf
        F = sf.F
        if self.verbose:
            print(f"got {sf} via reference spectrum!")
        return F

    def molecule(self, mol):
        """Return a calibrated molecule. dev-time only"""
        molecule = self.mdict.get(mol)
        if mol in self.mol_list:
            # if the molecule is in the matrix, give it sensitivity factors and
            # quantification coefficients for all of the masses in the matrix
            i = self.mol_list.index(mol)
            F_vec = self.F_mat[i, :]  # rows in F_mat correspond to molecules
            Q_vec = self.Q_mat[:, i]  # columns in Q_mat correspond to molecules
            F_i = dict(list(zip(self.mass_list, F_vec)))
            Q_i = dict(list(zip(self.mass_list, Q_vec)))
            molecule.F = F_i
            molecule.Q = Q_i
        else:
            mass = molecule.get_primary()
            F = self.get_F(mol, mass)
            molecule.F = {mass: F}
        return molecule

    @property
    def Q_mat(self):
        """The quantification matrix. Flux per signal in [mol/C]. Calculated once."""
        if self._Q_mat is not None:  # the comparison to None is necessary here.
            return self._Q_mat
        return self.calc_Q_mat()

    def calc_Q_mat(self):
        """Make the quantification matrix that can convert signal to flux.

        Calculates the quantification matrix, Q_mat. Columns in Q_mat correspond to the
        masses in self.mass_list, and rows in Q_mat correspond to the molecules in
        self.mol_list. Entries have units [mol/C].  Multiplying a vector with the
        signals for each of the masses in self.mass_list (i.e., signal) by Q_mat results
        in a vector with the flux of each molecule in self.mol_list (i.e., n_dot_vec).
        Q_mat is, beautifully, just the inverse of self.F_mat.

        Returns: 2-D matrix of floats: self.Q_mat = self.F_mat^-1.
        """
        F_mat = self.F_mat
        if self.N_mass == self.N_mol:
            try:
                Q_mat = np.linalg.inv(F_mat)  # <-- ! THE LINEAR ALGEBRA MAGIC IS HERE !
            except np.linalg.LinAlgError:
                raise SensitivityMatrixError(
                    "cannot take inverse of this square sensitivity matrix:\n"
                    + self.prints_F_mat()
                    + f"\nMade from mol_list={self.mol_list} "
                    f"and mass_list={self.mass_list}"
                )
        elif self.N_mass > self.N_mol:
            try:
                Q_mat = np.linalg.pinv(F_mat)  # <--- ! ...AND EXTRA MAGIC IS HERE !
            except np.linalg.LinAlgError:
                raise SensitivityMatrixError(
                    "cannot take inverse of this non-square sensitivity matrix:\n"
                    + self.prints_F_mat()
                    + f"\nMade from mol_list={self.mol_list} "
                    f"and mass_list={self.mass_list}"
                )
        else:
            raise ValueError(
                "Can't make quantification matrix for "
                f"{self.N_mol} mols with {self.N_mass} masses: "
                f"\nmass_list = {self.mass_list}"
                f"\nmol_list = {self.mol_list}"
            )
        self._Q_mat = Q_mat
        return Q_mat

    def print_Q_mat(self):
        """Print the active quantification matrix in decent-looking format."""
        print(np.round(self.Q_mat, decimals=2))

    def calc_signal(self, n_dot):
        """Calculate and return the signal given flux by dot'ing with self.F_mat

        Args:
            n_dot (dict): flux in [mol/s] signals for molecules in self.mol_list
                n_dot does not need to be complete - any missing molecules will be
                assumed to have zero flux.

        Returns:
            signal (dict): the predicted signal in [A] for each mass in self.mass_list
        """
        n_dot_vec = np.array([n_dot.get(mol, 0) for mol in self.mol_list])
        signal_vec = self.F_mat.dot(n_dot_vec)
        signal = dict(zip(self.mass_list, signal_vec))
        return signal

    def calc_n_dot(self, signals):
        """Calculate and return the flux given signal by dot'ing with self.Q_mat

        Args:
            signals (dict or SignalDict): MID or advanced MID signals for at least each
                mass in self.mass_list, in [A].

        Returns:
            n_dot (dict): the flux in [mol/s] of each molecule in self.mol_list
        """
        S_vec = np.array([signals[mass] for mass in self.mass_list])
        n_dot_vec = self.Q_mat.dot(S_vec)
        n_dot = dict(zip(self.mol_list, n_dot_vec))
        return n_dot


def STANDARD_T_OF_M(M):  # should maybe move to the constants module?
    """The standard transmission function, a function of m/z in atomic units"""
    return M**STANDARD_TRANSMISSION_EXPONENT


class SensitivityFit:
    """Class for describing and using a trend in measured sensitivity factors

    SensitivityFit has two distinct states. Fitted and not fitted. The property
    fitted is True/False correspondingly. SensitivityFit switches from not fitted to
    fitted first time that the method fit() is called.
    """

    def __init__(
        self,
        sensitivity_list,
        *,
        setting=None,
        alpha=None,
        f_fun=None,
        T_of_M=None,
        beta=None,
        E_ion=STANDARD_IONIZATION_ENERGY,
    ):
        """Initiate a SensitivityFit with a SensitivityList and optional starting params

        Typically SensitivityFit will be initated just with sl (the SensitivityList) and
        then fit with the fit() method. But, optionally, the fit can be given with the
        other __init__ parameters. This is the case if loaded by a calibration, in which
        case normally it will be initialized alpha and beta.

        The predicted value of F, under the present model, is:
            F = alpha * f
        where
            f = k * sigma^i * P^i(fragment) * T_of_M(M)
        where
            T_of_M(M) = M**beta
        As described in .../Industrial R&D/Quantification/Reports/MS_Theory_v1.0.pdf

        As of now, alpha and beta fully describe the fit. alpha is the ratio of F to
        f, and beta is the exponent in the transmission function used to calculate f.
        The option also exists to specify T_of_M directly, in anticipation of a
        future when T_of_M is not just M raised to an exponent, but something more
        sophisticated or an interpolation from a measured curve.
        k is just a normalization constant for which there's no reason to specify.

        Args:
            sensitivity_list (SensitivityList): the wrapped SensitivityFactors to fit
            alpha (float): Optional. The fitted value of F^O2_M32
            f_fun (function): Optional. A function of mol and mass returning the
                predicted relative sensitivity factor, f^{mol}_{mass}
            beta (float): Optional. The transmission-function exponent
        """
        if setting:
            sensitivity_list = sensitivity_list.filter(setting=setting)
        self.sl = sensitivity_list
        self.alpha = alpha
        self._f_fun = f_fun
        self._k = None  # the a normalization constant for f
        self._T_of_M = T_of_M
        self.beta = beta
        self.E_ion = E_ion
        self.mdict = MoleculeDict()
        self.__vecs = None  # vectors of F, f, and M, used internally
        self._alpha_0 = None  # The initial guess at the transmission function
        self._f_fun = f_fun

    def __eq__(self, other):
        """Return whether this SensitivityFit is equal to another"""
        return isclose(self.alpha, other.alpha) and isclose(self.beta, other.beta)

    def as_dict(self):
        """Dictionary representation of self"""
        self_as_dict = dict(
            # Note that saving self.sl combined the original sensitivity list with setting
            sensitivity_list=self.sl.as_dict(),
            alpha=self.alpha,
            beta=self.beta,
            # k=self.k, E_ion=self.E_ion
        )
        return self_as_dict

    def __repr__(self):
        """The SensitivityFit is described by its SensitivityFactors, now with f."""
        self_as_str = f"{self.__class__}(["
        for sf in self.sl:
            self_as_str += f"\n\t{sf},"
        self_as_str += "\n])"
        return self_as_str

    def reset(self, **kwargs):
        """Shortcut to dictate from outside the fit parameters, e.g. alpha and beta"""
        if "alpha" in kwargs:
            self.alpha = kwargs.pop("alpha")
        if "beta" in kwargs:
            self.beta = kwargs.pop("beta")
        for key, value in kwargs.items():
            print(f"setting fit.{key} to {value}")
            setattr(self, key, value)

    @property
    def _vecs(self):
        """List of [F_vec, f_vec, M_vec], vectors fo those properties in sf_list."""
        """Generate vectors of F, f, and M, and the initial guess at alpha"""
        F_vec = np.array([])
        f_vec = np.array([])
        M_vec = np.array([])  # for fitting beta

        for sf in self.sl:
            F = sf.F
            f = self.f_fun(sf.mol, sf.mass)
            M = sf.M

            F_vec = np.append(F_vec, F)
            f_vec = np.append(f_vec, f)
            M_vec = np.append(M_vec, M)

        I_sort = np.argsort(f_vec)

        f_vec, F_vec, M_vec = f_vec[I_sort], F_vec[I_sort], M_vec[I_sort]
        return f_vec, F_vec, M_vec

    @property
    def alpha_0(self):
        """Initial guess at alpha for fitting FIXME: like _vecs, abuse of @property?"""
        if not self._alpha_0:
            F_vec = np.array([])
            for sf in self.sl:
                if sf.mol == REFERENCE_MOLECULE and sf.mass == REFERENCE_MASS:
                    alpha_0 = sf.F
                    break
                F_vec = np.append(F_vec, sf.F)
            else:
                alpha_0 = np.mean(F_vec)
            self._alpha_0 = alpha_0
        return self._alpha_0

    @property
    def T_of_M(self):
        """The transmission function, a function of m/z in atomic units"""
        if not self._T_of_M:
            if self.beta:
                self.make_T_of_M()
        return self._T_of_M

    @property
    def fitted(self):
        """Whether the SensitivityFit has calculated its fitting function"""
        return self.alpha is not None and self.T_of_M is not None

    def make_T_of_M(self):
        """Prepare the transmission function from the transmission function exponent"""

        def T_of_M(M):
            return M**self.beta

        self._T_of_M = T_of_M

    def fit_beta(self):
        """Find the value of beta minimizing the rms error of F vs f

        The way it works is it actually fits a delta_beta which fixes the f's based
        on the existing (old or default) fit in order to give the best new fit. It then
        adds this delta_beta to self.beta

        It also updates self.k in the end.
        """
        f_vec, F_vec, M_vec = self._vecs  # this calls f_fun() on each (mol, mass)
        alpha_0 = self.alpha_0

        def square_error(params):
            """return the error of F-vs-f given alpha and a change in beta"""
            alpha_i, delta_beta_i = params
            # pred. sensitivity factors based on fs from outer scope and fit params:
            Fs_pred = alpha_i * f_vec * (M_vec**delta_beta_i)
            # compare with Fs from outer scope:
            rel_err_vec = 2 * (F_vec - Fs_pred) / (F_vec + Fs_pred)
            sq_rel_err = rel_err_vec.dot(rel_err_vec)
            return sq_rel_err

        delta_beta_0 = 0
        res = minimize(
            square_error,
            np.array([alpha_0, delta_beta_0]),
        )
        alpha_wrong, delta_beta = res.x
        # print(f"{self}.fit_beta() gives the following res: \n {res}")
        self.beta = self.beta + delta_beta
        self._T_of_M = None
        self._calc_k()
        return self.beta

    def _calc_k(self):
        """Calculate k such that f(REFERENCE_MOLECULE, REFERENCE_MASS) = 1"""
        T_of_M = self.T_of_M
        molecule = self.mdict.get(REFERENCE_MOLECULE)
        spectrum = molecule.calc_norm_spectrum()
        sigma = molecule.calc_sigma(E_ion=self.E_ion)
        M = float(REFERENCE_MASS[1:])
        f_ref_without_k = sigma * spectrum[REFERENCE_MASS] * T_of_M(M)
        k = 1 / f_ref_without_k
        self._k = k

    def f_fun(self, mol, mass):
        """Return the predicted relative sensitivity factor for mol at mass"""
        if not self.fitted:
            self.fit()
        k = self.k
        T_of_M = self.T_of_M
        molecule = self.mdict.get(mol)
        spectrum = molecule.calc_norm_spectrum()
        pure_mass = mass_to_pure_mass(mass)
        if pure_mass not in spectrum:
            # print(
            #     f"Predicting 0 sensitivity for {mol} at {mass} "
            #     f"because {mass} is not in {mol}'s spectrum"
            # )
            return 0
        sigma = molecule.calc_sigma(E_ion=self.E_ion)
        M = mass_to_M(mass)  # turns a string e.g. "M44" to a float eg 44.0
        f = k * sigma * spectrum[pure_mass] * T_of_M(M)
        return f

    @property
    def k(self):
        """Normalization constant setting f(REFERENCE_MOLECULE, REFERENCE_MASS) = 1"""
        if not self._k:
            self._calc_k()
        return self._k

    def update_fs(self):
        """Go through the sensitivity factors and update their f with self.f_fun

        This copies each of the SensitivityFactors, but it is the same SensitivityList
        """
        for i, sf in enumerate(self.sl):
            new_sf = sf.copy()
            new_sf.f = self.f_fun(sf.mol, sf.mass)
            self.sl.sf_list[i] = sf

    def fit_alpha(self):
        """Return and set self.alpha which minimizes error of F^i_M = alpha * f^i_M"""

        f_vec, F_vec, M_vec = self._vecs
        alpha_0 = self.alpha_0

        def square_error(alpha_i):
            """return the relative square error of F-vs-f given just alpha"""
            # predict sensitivity factors based on fs from outer scope and fit alpha:
            F_vec_predicted = alpha_i * f_vec
            # compare with Fs from outer scope:
            rel_err_vec = 2 * (F_vec - F_vec_predicted) / (F_vec + F_vec_predicted)
            sq_rel_err = rel_err_vec.dot(rel_err_vec)
            return sq_rel_err

        res = minimize(
            square_error,
            np.array([alpha_0]),
        )
        (alpha,) = res.x
        # print(f"{self}.fit_alpha() gives the following res: \n {res}")
        self.alpha = alpha
        return alpha

    def fit(self):
        """Fit both beta and alpha, doing all the preparation steps needed in between

        First, this function makes sure that there is a fits-guess of beta and alpha,
        then fits beta (see SensitivitiyFit.fit_beta) and thus the predicted relative
        sensitivity factor f, and then fits alpha (see SensitivityFit.fit_alpha) and is
        thus ready to predict the absolute sensitivity factor of any mol, mass.

        The preparation steps (generating vectors of F, f, and M in _vecs; calculating
        the normalization constant k) are now called by fit_beta and fit_alpha.

        After fitting, this method updates the f values in each of the SensitivityFit's
        SensitivityFactors according to the fit.
        """
        if not self.beta:
            self.beta = STANDARD_TRANSMISSION_EXPONENT
        if not self.alpha:
            self.alpha = self.alpha_0
        self.fit_beta()
        self.fit_alpha()
        self.update_fs()

    def predict_F(self, mol, mass):
        """Predict absolute sensitivity factor for mol at mass as float in [C/mol]"""
        f = self.f_fun(mol, mass)
        F = self.alpha * f
        return F

    def predict_sf(self, mol, mass):
        """Predict sensitivity factor and return as a SensitivityFactor instance"""
        f = self.f_fun(mol, mass)
        F = self.alpha * f
        return SensitivityFactor(mol=mol, mass=mass, F=F, f=f, F_type="predicted")

    def plot_F_vs_f(
        self,
        ax="new",
        predict=None,
        labels=True,
        plot_fit=False,
    ):
        """Plot active (measured) vs predicted sensitivity factors

        This is a way to visualize and sanity-check the active calibration factors
        stored in the calibration. Each active sensitivity factor (F^i_M) is plotted
        against a corresponding predicted sensitivity factor (f^i_M), where i is the
        molecule and M is the m/z value. F has the units [C/mol] and f is relative
        to F^O2_M32 and thereby unitless.
        Each point is colored with the standard EC-MS color for m/z=M, shaped with a
        symbol representing how the calibration, and (if labels) labeled with F^i_M.

        Args:
            ax (matplotlib...Axes...): The matplotlib axis handle to plot on.
                By default a new axis is created.
            predict (dict): {mol:mass} where you want the predicted sensitivity for
                mol at mass shown on the plot.
            labels (bool): Whether to add labels. Manual labels in Inkscape look best.
            plot_fit (bool): Whether to plot the fit line which is used when F_i_M has
                to be predicted.
        """
        sl = self.sl
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel(
                "f / f$^{" + REFERENCE_MOLECULE + "}_{" + REFERENCE_MASS + "}$"
            )
            ax.set_ylabel("F / [C/mol]")
        if predict is not None:
            for mol, mass in predict.items():
                sl.append(self.predict_sf(mol, mass))
        f_max = 0
        for sf in sl:
            mol = sf.mol
            mass = sf.mass
            F_type = sf.F_type
            spec = CAL_TYPE_SPECS.get(F_type, {"marker": "+"})
            try:
                color = STANDARD_MOL_COLORS[mol]
            except KeyError:
                color = STANDARD_COLORS.get(mass_to_pure_mass(mass), None)
            spec.update(color=color)
            label = "F$^{" + mol + "}_{" + mass + "}$"
            spec.update(label=label)
            F_i_M = sf.F
            f_i_M = self.f_fun(mol, mass)
            f_max = max([f_max, f_i_M])
            ax.plot(f_i_M, F_i_M, **spec)
            if labels:
                ax.annotate(label, xy=[f_i_M, F_i_M])
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        if plot_fit:
            f_fit = np.array([0, f_max])
            F_fit = self.alpha * f_fit
            ax.plot(f_fit, F_fit, "k--")
        return ax
