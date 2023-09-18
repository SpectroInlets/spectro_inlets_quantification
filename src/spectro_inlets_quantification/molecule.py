# ruff: noqa: E741

# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Everything to do with physical properties of molecules.

Variables with abbreviated or non-descriptive names, e.g. physical quantities:

* T: Temperature in [K]
* p: pressure in [Pa]
* M: molar mass in [g/mol] **NOTE**: Not SI!
* D: diffusion constant in [m^2/s] in water at standard conditions unless otherwise noted
* Hcp: solubility constant (concentration/pressure) in [(mol/l)/bar] **NOTE**: Not SI!
* H_0: Hcp at standard temperature in [(mol/l)/bar] **NOTE**: Not SI!
* T_c: temperature-dependence of Hcp in [K]
* KH: volatility constant (pressure/concentration) in [Pa/(mol/m^3)]
* sigma: ionization cross section in [Ang^2] **NOTE**: Not SI!
* beta: (fitted) exponent in the transmission function [dimensionless]

The variable names are from .../Industrial R&D/Quantification/Reports/MS_Theory_v1.0

"""
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Optional, Union, cast

import numpy as np
import yaml

from attrs import Attribute, asdict, define, field

from .compatability import TypeAlias
from .config import Config
from .constants import (
    AVOGADRO_CONSTANT,
    GAS_CONSTANT,
    PKAS,
    STANDARD_COLORS,
    STANDARD_IONIZATION_ENERGY,
    STANDARD_MOL_COLORS,
    STANDARD_TEMPERATURE,
)
from .custom_types import MASS
from .medium import Medium
from .tools import Singleton, make_axis

if TYPE_CHECKING:
    from matplotlib import pyplot

CONFIG = Config()

MASS_TO_FLOAT: TypeAlias = Dict[str, float]
PATH_OR_STR: TypeAlias = Union[Path, str]
TRANSMISSION_AMPLIFICATION_FUNCTION: TypeAlias = Callable[[float], float]


class MoleculeDict(dict, metaclass=Singleton):  # type: ignore
    """This class implements the MoleculeDict.

    It is a mapping of molecule names to Molecule, which loads the Molecule, if it hasn't been
    already.

    .. warning::
        This class overrides only `get` meaning that any methods that will create a new
        dictionary, will create an ordinary dict and not a MoleculeDict

    """

    def __init__(self, *args: Any, **kwargs: Any):
        """Initialize this objects attributes."""
        super().__init__(*args, **kwargs)
        self.medium = Medium()  # this is likely THE place medium is defined.

    def get(self, key: str) -> "Molecule":  # type: ignore
        """Return ``Molecule -> Molecule.load(key)`` first time, then stored `Molecule` later."""
        if key not in self:
            self[key] = Molecule.load(key)
        return cast("Molecule", self[key])


# ---------  the class ------------- #
@define(slots=False)
class Molecule:
    """The Molecule class. Interfaces with physical and chemical data.

    Args:
        name (str): Name of the molecule. Equal to the file name (minus extension)
            in the molecule data directory. Some are formulas, e.g. "CH4", and
            others are full names, e.g. "propionic acid".
        real_name (str): The real name of the molecule
        formula (str): Formula indicating which atoms are in the molecule
        spectrum (dict): {M: I_M} where I_M where M is a mass string (e.g. "M32")
            and I_M of is the (relative) intensity at that mass in the molecule's
            active spectrum. Allowed to change, unlike spectrum_0.
        spectrum_0 (dict): {M: I_M} for the reference spectrum.
        molecule_diameter (float): Gas-phase molecule diameter in [m]
        dynamic_viscosity (float): Gas-phase dynamic viscosity in [Pa*s]
        density_RTP (float): Gas-phase density at standard conditions in [kg/m^3]
        D_gas_RTP (float): Gas-phase diffusion constant at std. conds. in [m^2/s]
        M (float): Molecular mass in [g/mol]
        D (float): Liquid-phase diffusion constant at std. conds. in [m^2/s]
        H_0 (float): Henry's-Law solubility (concentration/pressure) constant at
            standard T in [(mol/l)/bar]
        T_c (float): Constant of Henry's-Law temperature dependence in [K]
        sigma (dict): {E_ion: sigma(E_ion)} where E_ion is ionization energy in
            [eV] and sigma(E_ion) is the ionization cross section in [AA^2]
        primary (str): Default mass for calibration (e.g. "M32")
        beta (float): Exponent to the transmission-amplification function
        E_ion (int): Internal condition - ionization energy
        thermo (dict): Themochemistry data including standard enthalpy of formation
            in [kJ/mol] and standard entropy in [J/(mol*K)] for various phases.
        kH (float): Obsolete, will be removed in a future version
        T_of_M (function): Transmission-amplification function, a function of mass
             by which to weigh intensities in the spectra. If one is not given, but beta is, a
             default transmission-amplification function will be created with the formula:
             T_of_M(M) = M ** beta
        verbose (bool): Whether to print stuff to the terminal

    """

    name: str = field()
    real_name: Optional[str] = field(default=None)
    formula: Optional[str] = field(default=None)
    spectrum: Optional[MASS_TO_FLOAT] = field(default=None)
    spectrum_0: Optional[MASS_TO_FLOAT] = field(default=None)
    molecule_diameter: Optional[float] = field(default=None)
    dynamic_viscosity: Optional[float] = field(default=None)
    density_RTP: Optional[float] = field(default=None)
    D_gas_RTP: Optional[float] = field(default=None)
    M: Optional[float] = field(default=None)
    D: Optional[float] = field(default=None)
    H_0: Optional[float] = field(default=None)
    T_c: Optional[float] = field(default=None)
    sigma: Dict[int, float] = field(default=None)
    primary: Optional[str] = field(default=None)
    beta: Optional[float] = field(default=None)
    E_ion: Optional[int] = field(default=STANDARD_IONIZATION_ENERGY)
    thermo: Optional[Dict[str, Dict[str, float]]] = field(default=None)
    kH: Optional[float] = field(default=None)
    T_of_M: Optional[TRANSMISSION_AMPLIFICATION_FUNCTION] = field(
        default=None, metadata={"serialize": False}
    )
    verbose: Optional[bool] = field(default=False, metadata={"serialize": False})

    # Initialize properties, which are not passed into init as arg
    corr_spectrum: Optional[MASS_TO_FLOAT] = field(init=False, default=None)
    H: Optional[float] = field(init=False, default=None)
    n_dot_0: Optional[float] = field(init=False, default=None)
    medium: Medium = field(init=False, factory=Medium)

    def __attrs_post_init__(self) -> None:
        """Update initialized values after ``__init__``."""
        if self.spectrum is None and self.spectrum_0 is not None:
            self.spectrum = self.spectrum_0.copy()

        if self.T_of_M is None and self.beta is not None:

            def T_of_M(M: float) -> float:
                return cast(float, M**self.beta)

            self.T_of_M = T_of_M

    def as_dict(self) -> Dict[str, Any]:
        """Return a dictionary including everything needed to recreate self."""

        def should_serialize(attribute: Attribute, _: Any) -> bool:  # type: ignore
            """Filter attributes that should not be serialized.

            Filter function that makes sure that only attributes that can be init'ed and
            have not been marked as don't serialize are included in as_dict version
            """
            return attribute.init and attribute.metadata.get("serialize", True)

        return asdict(self, filter=should_serialize)

    @cached_property
    def norm_spectrum(self) -> MASS_TO_FLOAT:
        """Normalized version of self.spectrum."""
        return self.calc_norm_spectrum()

    def save(self, file_name: Optional[str] = None, mol_dir: Optional[PATH_OR_STR] = None) -> None:
        """Save the `as_dict` form of the molecule to a yaml file.

        Saves in `CONFIG.aux_data_directory / "molecules"` if it is set (this is the
        user's quant data library, as opposed to the included library). If it is not set,
        saves in `CONFIG.data_directory / "molecules"`.

        Args:
            file_name: Name of the yaml file, including the file extension ".yml"
            mol_dir: Path to directory to save molecule in. Defaults are as outlinied above.
        """
        if file_name:
            file_name_with_suffix = Path(file_name).with_suffix(".yml")
        else:
            file_name_with_suffix = Path(self.name + ".yml")

        path_to_yaml = CONFIG.get_save_destination(
            data_file_type="molecules",
            filepath=file_name_with_suffix,
            override_destination_dir=mol_dir,
        )
        self_as_dict = self.as_dict()
        with open(path_to_yaml, "w") as yaml_file:
            yaml.dump(self_as_dict, yaml_file, indent=4)

    @classmethod
    def load(
        cls, file_name: str, mol_dir: Optional[PATH_OR_STR] = None, **kwargs: Any
    ) -> "Molecule":
        """Loads a molecule object from a yaml file.

        Args:
            file_name: Name of the yaml file WITHOUT the ".yml" extension
            mol_dir: Path to directory to save molecule in, defaults to
                :attr:`Config.molecule_directory` in order
            kwargs: (other) key word arguments are fed to Molecule.__init__()

        Returns:
            Molecule: a Molecule object ready to inform your calculations!
        """
        file_name_with_suffix = Path(file_name).with_suffix(".yml")
        try:
            file_path = CONFIG.get_best_data_file(
                data_file_type="molecules",
                filepath=file_name_with_suffix,
                override_source_dir=mol_dir,
            )
        except ValueError as value_error:
            raise ValueError(
                f"Can't find a molecule named '{file_name}'. Please consider providing an "
                "`aux_data_directory` (which contains a 'molecules' folder) to "
                "`config.Config` or a ``mol_dir`` to this method, either of which contains "
                "the molecule file."
            ) from value_error

        with open(file_path) as yaml_file:
            self_as_dict = yaml.safe_load(yaml_file)

        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    def update(self, **kwargs: Any) -> None:
        """Set attributes given as kwargs, but update rather than replacing dicts."""
        for key, value in kwargs.items():
            if isinstance(value, dict):
                try:
                    getattr(self, key).update(value)
                except AttributeError:
                    setattr(self, key, value)
                continue
            setattr(self, key, value)

    @property
    def eta(self) -> Optional[float]:
        """Pseudonym for 'dynamic_viscosity' in [Pa*s]."""
        return self.dynamic_viscosity

    @property
    def s(self) -> Optional[float]:
        """Pseudonym for 'molecular_diameter' in [m]."""
        return self.molecule_diameter

    @property
    def m(self) -> float:
        """The molecule mass in [kg]."""
        M = self.M
        m = M * 1e-3 / AVOGADRO_CONSTANT  # from g/mol to kg
        return m

    def get_primary(self) -> MASS:
        """Return the default mass for quantification.

        Pre-defined as ``self.primary`` or max of spectrum.
        """
        if self.primary is not None:
            return self.primary
        if self.spectrum is not None:
            masses, values = zip(*list(self.spectrum.items()))
            index = np.argmax(values)
            return masses[index]

    def calc_sigma(self, E_ion: Optional[int] = None) -> float:
        """Return the ionization cross-section [A^2] given the ionization energy [eV].

        Checks if E_ion (e.g. 80) is in the keys of self.sigma (e.g. {70:1.5, 100:4.5}),
        returns self.sigma[E_ion] if it is and otherwise interpolates (e.g. getting 2.5)

        Args:
            E_ion (float): Ionization energy in [eV]

        Raises:
            TypeError: If interpolation fails
            AttributeError: If interpolation results in a negative value

        Returns:
            float: The ionization cross-section in [A^2]
        """
        sigma = self.sigma
        if sigma is None:
            raise AttributeError(f"Can't get E_ion={E_ion} from {self.name} with sigma={sigma}")
        if E_ion is None:
            E_ion = self.E_ion
        if E_ion in sigma:
            return sigma[E_ion]
        else:
            if self.verbose:
                print(f"getting sigma at E_ion = {E_ion} for {self.name} by interpolation")
            E, sig = zip(*list(sigma.items()))
            E, sig = np.array([float(E_n) for E_n in E]), np.array(sig)
            I_sort = np.argsort(E)
            E, sig = E[I_sort], sig[I_sort]
            if E_ion < min(E) or E_ion > max(E):
                if self.verbose:
                    print(
                        f"Warning! You want E_ion={E_ion} from {self.name} which has "
                        f"sigma={sigma}. Using data endpoint nearest requested value."
                    )
            try:
                sigma_single = cast(
                    float,
                    np.interp(
                        E_ion,
                        E,
                        sig,
                    ),
                )
            except TypeError:
                print(f"Can't get E_ion={E_ion} from {self.name} with sigma={sigma}")
                raise TypeError
            if sigma_single < 0:
                raise AttributeError(f"Can't get E_ion={E_ion} from {self.name} with sigma={sigma}")

            return sigma_single

    def calc_norm_spectrum(self) -> MASS_TO_FLOAT:
        """Return and set the normalized active spectrum as a dict {M:I_M}."""
        spectrum = self.spectrum
        total_intensity = sum(list(spectrum.values()))
        norm_spectrum = {}
        for mass, I in spectrum.items():
            norm_spectrum[mass] = I / total_intensity
        self.norm_spectrum = norm_spectrum
        return norm_spectrum

    def correct_spectrum(
        self,
        T_of_M: Optional[TRANSMISSION_AMPLIFICATION_FUNCTION] = None,
        beta: Optional[float] = None,
    ) -> None:
        """Set self.spectrum to a transmission-amplification-weighted spectrum.

        Args:
            T_of_M (function): The transmission-amplification function, defaults to self.T_of_M
                if not set
            beta (float): The exponent to the transmission-amplification function,
                used if T_of_M is not given. NOTE: This values will be set via `set_beta`,
                which generates a new ``self.T_of_M`` even if T_of_M is given as an argument.

        Returns:
            dict: {M:I_M} where I_M is intensity at mass M in the corrected spectrum
        """
        if beta is not None:
            self.set_beta(beta)
        if T_of_M is None:
            T_of_M = self.T_of_M
        else:
            self.T_of_M = T_of_M
        self.spectrum = self.calc_corr_spectrum(T_of_M=T_of_M)

    def set_beta(self, beta: float) -> None:
        """Set transmission-amplification function via its exponent."""
        self.beta = beta

        def T_of_M(M: float) -> float:
            return cast(float, M**beta)

        self.T_of_M = T_of_M

    def calc_corr_spectrum(
        self, T_of_M: Optional[TRANSMISSION_AMPLIFICATION_FUNCTION] = None
    ) -> MASS_TO_FLOAT:
        """Correct the spectrum by weighing by a transmission-amplification function.

        This function stores the corrected spectrum as self.corr_spectrum, doesn't touch
        self.spectrum.

        Args:
            T_of_M (function): Transmission-amplification function

        Returns:
            dict: {M: I_M} with I_M being the intensity at mass M in the normalized and
            corrected spectrum
        """
        if T_of_M is None:
            T_of_M = self.get_T_of_M()

        spectrum = {}
        for mass, I in self.spectrum_0.items():
            # ^ Important to use spectrum_0, to not accumulate beta's
            M = float(mass[1:])
            spectrum[mass] = I * T_of_M(M)
        total_intensity = sum(list(spectrum.values()))
        corr_spectrum = {}
        for mass, I in spectrum.items():
            corr_spectrum[mass] = I / total_intensity
        self.corr_spectrum = corr_spectrum
        return corr_spectrum

    def plot_spectrum(
        self,
        norm: Optional[bool] = False,
        T_of_M: Optional[TRANSMISSION_AMPLIFICATION_FUNCTION] = None,
        top: Optional[float] = 1,
        width: Optional[float] = 0.5,
        offset: Optional[float] = 0,
        ax: Union[str, "pyplot.Axes"] = "new",
        **kwargs: Any,
    ) -> "pyplot.Axes":  # pragma: no cover
        """Plots the molecule's reference mass spectrum.

        Args:
            norm (bool): Whether to normalize to the sum of the peaks
            T_of_M (function): Transmission-amplification function
            top (float): Height of the highest peak in the plot
            width (float): Width of the bars in the plot (arg to matplotlib.pyplot.bar)
            offset (float): How much to shift bars on m/z axis (useful for co-plotting)
            ax (str or matplotlib Axis): Axis to use. Default "new" creates a new one.
            kwargs (str): Key word arguments passed to (arg to matplotlib.pyplot.bar)

        Returns:
            matplotlib Axes: The axes on which the spectrum was plotted
        """
        if T_of_M is not None:
            print(f"plotting transmission-corrected spectrum for {self.name}")
            spectrum = self.calc_corr_spectrum(T_of_M=T_of_M)
        else:
            print(f"plotting un-corrected spectrum for {self.name}")
            spectrum = self.calc_norm_spectrum()
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("m/z")
            ax.set_ylabel("norm. intensity")
            ax.set_title(f"{self.name} reference spectrum")
        ax = cast("pyplot.Axes", ax)
        if norm:
            factor = 1 / sum(spectrum.values())
        else:
            factor = top / max(spectrum.values())
        if "color" not in kwargs:
            kwargs["color"] = STANDARD_MOL_COLORS.get(self.name, None)
        for mass, value in spectrum.items():
            M = float(mass[1:])
            ax.bar(M + offset, value * factor, width=width, **kwargs)
        return ax

    def get_color(self) -> str:
        """Return the molecule's color = the EC-MS standard color of its primary mass.

        Standard colors are stored in `constants.STANDARD_COLORS`

        """
        if self.name in STANDARD_MOL_COLORS:
            return STANDARD_MOL_COLORS[self.name]
        primary = self.get_primary()
        return STANDARD_COLORS.get(primary, "k")

    def get_T_of_M(self) -> TRANSMISSION_AMPLIFICATION_FUNCTION:
        """Return the active transmission-amplification function."""
        if hasattr(self, "T_of_M") and self.T_of_M is not None:
            return self.T_of_M
        elif hasattr(self, "beta") and self.beta is not None:

            def T_of_M(M: float) -> float:
                return cast(float, M**self.beta)

            return T_of_M
        else:
            raise AttributeError(
                f"Molecule {self.name} has no attr 'T_of_M' or 'beta', or both are None"
            )

    def calc_Hcp(self, T: Optional[float]) -> float:  # noqa: C901
        """Returns the solubility Henry's-Law constant in Sanders units: [(mol/l) / bar].

        Solubility is also called concentration/pressure (cp), thus the cp in Hcp.
        To get volatility (pc) in SI units: Take the reciprocal and multiply by 100, or use
        Molecule.calc_KH(T) instead.
        This function uses data copied over to the molecule file from Sanders'
        Henry's-Law compilation.

        Args:
            T (float): Temperature in [K], ``self.medium.T`` by default.

        Returns:
            float: Hcp in [(mol/l) / bar]
        """
        if T is None:
            T = self.medium.T
        # Example: 1.2E-03 * EXP( 1700.*(1./298.15-1./T) )
        T_c = self.T_c
        H_0 = self.H_0
        if T_c is None:
            if self.verbose:
                print(f"Warning!!! T_c={T_c} for mol={self.name}")
            try:
                dfH0 = self.thermo["dfH0"]
                d_vap_H = (dfH0["gas"] - dfH0["liquid"]) * 1e3
            except KeyError:
                if self.verbose:
                    print("Assuming no temperature dependence")
                T_c = 0
            else:
                if self.verbose:
                    print(r"Using ($\Delta_{vap}$H/R instead!) :D ")
                T_c = d_vap_H / GAS_CONSTANT
        if H_0 is None:
            if self.verbose:
                print(f"Warning!!! H_0={H_0} for mol={self.name}")
            try:
                dfH0 = self.thermo["dfH0"]
                S0 = self.thermo["S0"]
                d_vap_H = (dfH0["gas"] - dfH0["liquid"]) * 1e3  # in J/mol
                d_vap_S = S0["gas"] - S0["liquid"]  # in J/(mol*K)
            except KeyError:
                if self.verbose:
                    print("Assuming zero volatility")
                H_0 = 0
            else:
                if self.verbose:
                    print(r"Using ($\Delta_{vap}$G) and the molar density instead")
                d_vap_G = d_vap_H - STANDARD_TEMPERATURE * d_vap_S  # in J/mol
                if hasattr(self, "rho_l"):
                    rho_l: float = self.rho_l  # type: ignore
                else:
                    if self.verbose:
                        print("Assuming 1000 kg/m^3 (same density as water)")
                    rho_l = 1e3
                c_0 = rho_l * 1e3 / self.M  # g/m^3 / [g/mol] = mol/m^3
                H_0 = (
                    c_0 * np.exp(d_vap_G / (GAS_CONSTANT * STANDARD_TEMPERATURE)) * 1e-3
                )  # (mol/m^3)/bar * m^3/l
                # sign: more positive d_vap_G means it likes being liquid,

        H = H_0 * cast(float, np.exp(T_c * (1.0 / T - 1.0 / STANDARD_TEMPERATURE)))
        return H

    def calc_KH(self, T: Optional[float] = None) -> float:
        """Return the volatility Henry's-Law constant in SI units [Pa/[mol/m^3]].

        T will default to ``self.medium.T`` if not given

        """
        Hcp = self.calc_Hcp(T=T)  # in M/bar
        KH = 100 / Hcp  # in Pa/mM = ([Pa/bar]*[M/mM]) / [M/bar]
        # ^ where the unit converter is (100 [Pa/bar]*[M/mM]) = 1
        return KH

    def calc_H(
        self,
        n_dot_0: Optional[float],
        p: Optional[float] = None,
        T: Optional[float] = None,
    ) -> float:
        """Return the molecule's mass-transfer number in [m^3/s].

        Args:
            n_dot_0 (float): The total flux through the capillary of the chip in [mol/s]
            p (float): The pressure in [Pa], defaults to ``self.medium.p`` if not given
            T (float): The temperature in [K], defaults to ``self.medium.T`` if not given

        Returns:
            float: Mass-transfer number, i.e. ratio of flux to concentration in [m^3/s]
        """
        if n_dot_0 is None:
            n_dot_0 = self.n_dot_0
        else:
            self.n_dot_0 = n_dot_0
        KH = self.calc_KH(T=T)
        p = p or self.medium.p
        H = KH * n_dot_0 / p
        # [m^3/s] = [Pa/(mol/m^3)] * [mol/s] / [Pa]
        self.H = H
        return H

    def calc_p_vap(self, T: Optional[float] = None) -> float:
        """Return the vapor pressure of the molecule in [Pa] given temperature in [K].

        T defaults to ``self.medium.T`` if not given

        """
        if T is None:
            T = self.medium.T
        try:
            dfH0 = self.thermo["dfH0"]
            S0 = self.thermo["S0"]
            dH = (dfH0["gas"] - dfH0["liquid"]) * 1e3
            dS = S0["gas"] - S0["liquid"]
        except KeyError:
            print(
                f"{self.name} does not have the "
                "thermochem data needed to calculate p_vap! Returning 0."
            )
            return 0

        p0 = 1e5  # [Pa]

        p_vap = p0 * cast(float, np.exp(-dH / (GAS_CONSTANT * T) + dS / GAS_CONSTANT))

        return p_vap

    @property
    def pKa(self) -> float:
        """Return the pKa above or below which (see pKa_description) self is volatile."""
        pKa, description = PKAS.get(self.name, (None, None))
        if self.verbose:
            if pKa:
                print(f"{self.name} '{description}' pKa={pKa}")
            else:
                print(f"{self.name} does not have a recorded pKa. PKAs = {PKAS}")
        return pKa

    @property
    def pKa_description(self) -> str:
        """Return str specifying whether the volatile compound is at pH below or above pKa."""
        pKa, description = PKAS.get(self.name, (None, None))
        return description

    def calc_volatile_portion(self, pH: float) -> float:
        """Return the fraction in the volatile form of mol as a function of pH."""
        pKa, description = PKAS.get(self.name, (None, None))
        HA_to_A_ratio = cast(float, np.power(10, pKa - pH))
        portion_HA = HA_to_A_ratio / (1 + HA_to_A_ratio)
        if description == "volatile below":
            return portion_HA
        if description == "volatile above":
            return 1 - portion_HA
        raise NotImplementedError(f"Don't know how {self.name} can be '{description}' pKa={pKa}")
