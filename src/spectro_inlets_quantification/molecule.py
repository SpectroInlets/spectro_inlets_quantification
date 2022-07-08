# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Everything to do with physical properties of molecules

Variables with abbreviated or non-descriptive names, e.g. physical quantities:

    T: Temperature in [K]
    p: pressure in [Pa]
    M: molar mass in [g/mol]                                            # not SI!
    D: diffusion constant in [m^2/s]
        in water at standard conditions unless otherwise noted
    Hcp: solubility constant (concentration/pressure) in [(mol/l)/bar]  # not SI!
    H_0: Hcp at standard temperature in [(mol/l)/bar]                   # not SI!
    T_c: temperature-dependence of Hcp in [K]
    KH: volatility constant (pressure/concentration) in [Pa/(mol/m^3)]
    sigma: ionization cross section in [Ang^2]                          # not SI!
    beta: (fitted) exponent in the transmission function [dimensionless]

The variable names are from .../Industrial R&D/Quantification/Reports/MS_Theory_v1.0

"""
import json
from pathlib import Path
import numpy as np
from .constants import (
    STANDARD_COLORS,
    STANDARD_TEMPERATURE,
    GAS_CONSTANT,
    STANDARD_IONIZATION_ENERGY,
    AVOGADRO_CONSTANT,
    PKAS,
    STANDARD_MOL_COLORS,
)
from .tools import make_axis, Singleton
from .medium import Medium
from .config import Config


CONFIG = Config()


class MoleculeDict(dict, metaclass=Singleton):
    def __init__(self, *args, medium=None, **kwargs):
        super().__init__(*args, **kwargs)
        if not medium:
            medium = Medium()  # this is likely THE place medium is defined.
        self.medium = medium

    def get(self, key):  # Using this, a Molecule is loaded if not already there
        """Return Molecule: Molecule.load(key) first time, the stored Molecule later."""
        if key not in self:
            self[key] = Molecule.load(key, medium=self.medium)
        return self[key]


# ---------  the class ------------- #
class Molecule:
    """The Molecule class. Interfaces with physical and chemical data."""

    def __init__(
        self,
        *,
        name=None,
        formula=None,
        spectrum=None,
        spectrum_0=None,
        molecule_diameter=None,
        dynamic_viscosity=None,
        density_RTP=None,
        D_gas_RTP=None,
        M=None,
        D=None,
        H_0=None,
        T_c=None,
        sigma=None,
        T=None,
        p=None,
        primary=None,
        T_of_M=None,
        beta=None,
        E_ion=STANDARD_IONIZATION_ENERGY,
        thermo=None,
        medium=None,
        verbose=False,
        **kwargs,
    ):
        """Create a Molecule object given its properties

        Args:
            name (str): Name of the molecule. Equal to the file name (minus extension)
                in the molecule data directory. Some are formulas, e.g. "CH4", and
                others are full names, e.g. "propionic acid". TODO: consistency on that?
            formula (str): Formula indicating which atoms are in the molecule
            spectrum (dict): {M: I_M} where I_M where M is a mass string (e.g. "M32")
                and I_M of is the (relative) intensity at that mass in the molecule's
                active spectrum. Allowed to change, unlike spectrum_0.
            spectrum_0 (dict): {M: I_M} for the reference spectrum.
            molecule_diameter (float): gas-phase molecule diameter in [m]
            dynamic_viscosity (float): gas-phase dynamic viscosity in [Pa*s]
            density_RTP (float): gas-phase density at standard conditions in [kg/m^3]
            D_gas_RTP (float): gas-phase diffusion constant at std. conds. in [m^2/s]
            M (float): molecular mass in [g/mol] TODO: Change this from common to SI?
            D (float): liquid-phase diffusion constant at std. conds. in [m^2/s]
            H_0 (float): Henry's-Law solubility (concentration/pressure) constant at
                standard T in [(mol/l)/bar] TODO: Change this from common to SI units?
            T_c (float): Constant of Henry's-Law temperature dependence in [K]
            sigma (dict): {E_ion: sigma(E_ion)} where E_ion is ionization energy in
                [eV] and sigma(E_ion) is the ionization cross section in [AA^2]
            T (float): external condition - temperature in [K]
            p (float): external condition - pressure in [Pa]
            primary (str): default mass for calibration (e.g. "M32")
            T_of_M (function): transmission-amplification function, a function of mass
                by which to weigh intensities in the spectrum
            beta (float): exponent to the transmission-amplification function
            E_ion (float): internal condition - ionization energy
            thermo (dict): themochemistry data including standard enthalpy of formation
                in [kJ/mol] and standard entropy in [J/(mol*K)] for various phases.
            medium (Medium): stores the p and T
            verbose (bool): whether to print stuff to the terminal
            kwargs: additional kwargs are set as attributes.
        """
        # basic stuff
        self.name = name
        self.formula = formula
        self.M = M

        # gas transport properties
        self.molecule_diameter = molecule_diameter
        self.dynamic_viscosity = dynamic_viscosity
        self.density_RTP = density_RTP
        self.D_gas_RTP = D_gas_RTP

        # liquid phase transport properties
        self.D = D
        self.H_0 = H_0
        self.T_c = T_c

        # thermochemistry properties
        self.thermo = thermo

        # ionization-fragmentation properties:
        self.sigma = sigma
        self.spectrum = spectrum  # spectrum may be corrected
        self.norm_spectrum = self.calc_norm_spectrum()  # normalizes self.spectrum
        self.corr_spectrum = None  # to be set later
        if spectrum_0 is None:
            spectrum_0 = spectrum.copy()
        self.spectrum_0 = spectrum_0

        # mass spec options (can be in the file)
        self.primary = primary

        # mass spec options (shouldn't be in the file)
        self.T_of_M = T_of_M
        self.beta = beta
        self.E_ion = E_ion

        # mass transport properties that involve knowledge of outside world
        self.H = None  # mass transfer number [m^3/s], possibly set/used later
        self.n_dot_0 = None  # total capillary flux [mol/s], possibly set/used later

        if not medium:
            medium = Medium(p=p, T=T)
        self.medium = medium

        self.verbose = verbose

        # whatever I missed
        for key, value in kwargs.items():
            if not hasattr(self, key):
                setattr(self, key, value)
            else:
                raise RuntimeError(
                    f"you tried to set {key} while initializing a Molecule, but it "
                    f"already has self.{key} = {getattr(self, key)}"
                )

    def as_dict(self):
        """Return a dictionary including everything needed to recreate self."""
        self_as_dict = {}

        # basic stuff:
        self_as_dict.update(name=self.name, formula=self.formula, M=self.M)

        # gas transport properties:
        self_as_dict.update(
            molecule_diameter=self.molecule_diameter,
            dynamic_viscosity=self.dynamic_viscosity,
            density_RTP=self.density_RTP,
            D_gas_RTP=self.D_gas_RTP,
        )

        # liquid phase transport properties:
        self_as_dict.update(D=self.D, H_0=self.H_0, T_c=self.T_c)

        # thermochemistry properties:
        self_as_dict.update(thermo=self.thermo)

        # ionization-fragmentation properties:
        self_as_dict.update(sigma=self.sigma, spectrum=self.spectrum)

        # mass spec options (can be in the file):
        self_as_dict.update(primary=self.primary)

        return self_as_dict

    def save(self, mol_dir=None, file_name=None):
        """save the self.as_dict() form of the molecule to a .json file

        Args:
            file_name: name of the .json file. filename.endswith(".json")
            mol_dir: path to directory to save molecule in, defaults to
                :attr:`Config.molecule_directory`
        """
        mol_dir = mol_dir or CONFIG.molecule_directory
        if file_name is None:
            file_name = self.name + ".json"
        path_to_json = Path(mol_dir) / file_name
        self_as_dict = self.as_dict()
        with open(path_to_json, "w") as json_file:
            json.dump(self_as_dict, json_file, indent=4)

    @classmethod
    def load(cls, file_name, mol_dir=None, **kwargs):
        """loads a chip object from a .json file

        Args:
            file_name: name of the .json file. filename.endswith(".json")
            mol_dir: path to directory to save molecule in, defaults to
                :attr:`Config.molecule_directory`
            kwargs: (other) key word arguments are fed to Molecule.__init__()

        Returns:
            Molecule: a Molecule object ready to inform your calculations!
        """
        mol_dir = mol_dir or CONFIG.molecule_directory
        path_to_json = Path(mol_dir) / (file_name + ".json")
        with open(path_to_json) as json_file:
            self_as_dict = json.load(json_file)
        self_as_dict.update(kwargs)
        # Unfortunately, saving and loading a dict with integer keys to json
        #   turns the keys into strings. This fixes that:
        try:
            self_as_dict["sigma"] = dict(
                [(int(key), value) for key, value in self_as_dict["sigma"].items()]
            )
        except AttributeError:
            print(f"Warning!!! {file_name} has sigma={self_as_dict['sigma']}.")
        return cls(**self_as_dict)

    def update(self, **kwargs):
        """set attributes given as kwargs, but update rather than replacing dicts"""
        for key, value in kwargs.items():
            if isinstance(value, dict):
                try:
                    getattr(self, key).update(value)
                except AttributeError:
                    setattr(self, key, value)
                continue
            setattr(self, key, value)

    @property
    def T(self):
        return self.medium.T

    @T.setter
    def T(self, T):
        self.medium.T = T

    @property
    def p(self):
        return self.medium.p

    @p.setter
    def p(self, p):
        self.medium.p = p

    @property
    def eta(self):
        """Pseudonymn for dynamic_viscosity in [Pa*s]"""
        return self.dynamic_viscosity

    @property
    def s(self):
        """Pseudonym for molecular_diameter in [m]"""
        return self.molecule_diameter

    @property
    def m(self):
        """The molecule mass in [kg]"""
        M = self.M
        m = M * 1e-3 / AVOGADRO_CONSTANT  # from g/mol to kg
        return m

    def get_primary(self):
        """Return the default mass for quantification: pre-defined or max of spectrum"""
        if self.primary is not None:
            return self.primary
        if self.spectrum is not None:
            masses, values = zip(*list(self.spectrum.items()))
            index = np.argmax(values)
            return masses[index]

    def calc_sigma(self, E_ion=None):
        """Return the ionization cross-section [A^2] given the ionization energy [eV]

        TODO: there's some naming and type confusion here. Molecule.calc_sigma returns a
            float, but Molecule.sigma is a dict with numbers as keys (integer in
            practice, but float in principle). I think Molecule.sigma should change to a
            list of tuples and be renamed. Docstring to be fixed after decision on this.

        Checks if E_ion (e.g. 80) is in the keys of self.sigma (e.g. {70:1.5, 100:4.5}),
        returns self.sigma[E_ion] if it is and otherwise interpolates (e.g. getting 2.5)

        Args:
            E_ion (float): ionization energy in [eV]

        Raises:
            TypeError if interpolation fails
            AttributeError if interpolation results in a negative value

        Returns:
            float: The ionization cross-section in [A^2]
        """
        sigma = self.sigma
        if sigma is None:
            raise AttributeError(
                f"Can't get E_ion={E_ion} from {self.name} with sigma={sigma}"
            )
        if E_ion is None:
            E_ion = self.E_ion
        if E_ion in sigma:
            return sigma[E_ion]
        else:
            if self.verbose:
                print(
                    f"getting sigma at E_ion = {E_ion} for {self.name} by interpolation"
                )
            E, sig = zip(*list(sigma.items()))
            E, sig = np.array([float(E_n) for E_n in E]), np.array(sig)
            I_sort = np.argsort(E)
            E, sig = E[I_sort], sig[I_sort]
            if E_ion < min(E) or E_ion > max(E):
                if self.verbose:
                    print(
                        f"Warning! You want E_ion={E_ion} from {self.name} which has "
                        f"sigma={sigma}. Doing my best."
                    )
            try:
                sigma = np.interp(
                    E_ion,
                    E,
                    sig,
                )  # left=-1, right=-1  #
            except TypeError:
                print(f"Can't get E_ion={E_ion} from {self.name} with sigma={sigma}")
                raise TypeError
            if sigma < 0:
                raise AttributeError(
                    f"Can't get E_ion={E_ion} from {self.name} with sigma={self.sigma}"
                )

            return sigma

    def calc_norm_spectrum(self):
        """Return and set the normalized active spectrum as a dict {M:I_M}"""
        spectrum = self.spectrum
        total_intensity = sum(list(spectrum.values()))
        norm_spectrum = {}
        for mass, I in spectrum.items():
            norm_spectrum[mass] = I / total_intensity
        self.norm_spectrum = norm_spectrum
        return norm_spectrum

    def correct_spectrum(self, T_of_M=None, beta=None):
        """Set self.spectrum to a transmission-amplification-weighted spectrum

        Args:
            T_of_M (function): The transmission-amplification function
            beta (float): The exponent to the transmission-amplification function,
                used if T_of_M is not given.

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

    def set_beta(self, beta):
        """Set transmission-amplification function via its exponent"""
        self.beta = beta

        def T_of_M(M):
            return M**beta

        self.T_of_M = T_of_M

    def calc_corr_spectrum(self, T_of_M=None):
        """Correct the spectrum by weighing by a transmission-amplification function

        For an explanation and usage example, see slides 28-29 of:
        .../Quantification/Reports/20D14_Quantification_Halfway_Report.pptx

        This function stores the corrected spectrum as self.corr_spectrum, doesn't touch
        self.spectrum.

        Args:
            T_of_M (function): transmission-amplification function

        Returns:
            dict: {M: I_M} with I_M being the intensity at mass M in the normalized and
                corrected spectrum.
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
        self, norm=False, T_of_M=None, top=1, width=0.5, offset=0, ax="new", **kwargs
    ):
        """Plots the molecule's reference mass spectrum

        Args:
            norm (bool): whether to normalize to the sum of the peaks
            T_of_M (function): transmission-amplification function
            top (float): height of the highest peak in the plot
            width (float): width of the bars in the plot (arg to matplotlib.pyplot.bar)
            offset (float): How much to shift bars on m/z axis (useful for co-plotting)
            ax (str or matplotlib Axis): axis to use. Default "new" creates a new one.
            kwargs (str): key word arguments passed to (arg to matplotlib.pyplot.bar)

        Returns:
            matplotlib Axes: the axes on which the spectrum was plotted
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

    def get_color(self):
        """Return the molecule's color = the EC-MS standard color of its primary mass"""
        if self.name in STANDARD_MOL_COLORS:
            return STANDARD_MOL_COLORS[self.name]
        primary = self.get_primary()
        return STANDARD_COLORS.get(primary, "k")

    def get_T_of_M(self):
        """Return the active transmission-amplification function"""
        if hasattr(self, "T_of_M") and self.T_of_M is not None:
            return self.T_of_M
        elif hasattr(self, "beta") and self.beta is not None:

            def T_of_M(M):
                return M**self.beta

            return T_of_M
        else:
            raise AttributeError(
                f"Molecule {self.name} has no attr 'T_of_M' or 'beta', or both are None"
            )

    def calc_Hcp(self, T):
        """Returns the solubility Henry's-Law constant in Sanders units: [(mol/l) / bar]

        Solubility is also called concentration/pressure (cp), thus the cp in Hcp.
        To get volatility (pc) SI units: Take the reciprocal and multiply by 100, or use
        Molecule.calc_KH(T) instead.
        This function uses data copied over to the molecule file from Sanders'
        Henry's-Law compilation.

        Args:
            T (float): Temperature in [K]. By default uses self.T

        Returns:
            float: Hcp in [(mol/l) / bar]
        """
        if T is None:
            T = self.T
        # Example: 1.2E-03 * EXP(   1700.*(1./298.15-1./T))
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
                    print("assuming no temperature dependence")
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
                    print("assuming zero volatility")
                H_0 = 0
            else:
                if self.verbose:
                    print(r"using ($\Delta_{vap}$G) and the molar density instead")
                d_vap_G = d_vap_H - STANDARD_TEMPERATURE * d_vap_S  # in J/mol
                if hasattr(self, "rho_l"):
                    rho_l = self.rho_l
                else:
                    if self.verbose:
                        print("assuming 1000 kg/m^3 (same density as water)")
                    rho_l = 1e3
                c_0 = rho_l * 1e3 / self.M  # g/m^3 / [g/mol] = mol/m^3
                H_0 = (
                    c_0 * np.exp(d_vap_G / (GAS_CONSTANT * STANDARD_TEMPERATURE)) * 1e-3
                )  # (mol/m^3)/bar * m^3/l
                # sign: more positive d_vap_G means it likes being liquid,

        H = H_0 * np.exp(T_c * (1.0 / T - 1.0 / STANDARD_TEMPERATURE))
        return H

    def calc_KH(self, T=None):
        """Return the volatility Henry's-Law constant at in SI units [Pa/[mol/m^3]]"""
        Hcp = self.calc_Hcp(T=T)  # in M/bar
        KH = 100 / Hcp  # in Pa/mM = ([Pa/bar]*[M/mM]) / [M/bar]
        # ^ where the unit converter is (100 [Pa/bar]*[M/mM]) = 1
        return KH

    def calc_H(self, n_dot_0, p=None, T=None):
        """Return the molecule's mass-transfer number in [m^3/s]

        Args:
            n_dot_0 (float): the total flux through the capillary of the chip in [mol/s]
            p (float): The pressure in [Pa]
            T (float): The temperature in [K]

        Returns:
            float: mass-transfer number, i.e. ratio of flux to concentration in [m^3/s]
        """
        if n_dot_0 is None:
            n_dot_0 = self.n_dot_0
        else:
            self.n_dot_0 = n_dot_0
        KH = self.calc_KH(T=T)
        p = p or self.p
        H = KH * n_dot_0 / p
        # [m^3/s] = [Pa/(mol/m^3)] * [mol/s] / [Pa]
        self.H = H
        return H

    def calc_p_vap(self, T=None):
        """Return the vapor pressure of the molecule in [Pa] given temperature in [K]"""
        if T is None:
            T = self.T
        try:
            dfH0 = self.thermo["dfH0"]
            S0 = self.thermo["S0"]
            dH = (dfH0["gas"] - dfH0["liquid"]) * 1e3
            dS = S0["gas"] - S0["liquid"]
        except KeyError:
            print(
                f"{self.name} does not have the"
                + "thermochem data needed to calculate p_vap! Returning 0."
            )
            return 0

        p0 = 1e5  # [Pa]

        p_vap = p0 * np.exp(-dH / (GAS_CONSTANT * T) + dS / GAS_CONSTANT)

        return p_vap

    @property
    def pKa(self):
        """Return the pKa above or below which (see pKa_description) self is volatile"""
        pKa, description = PKAS.get(self.name, (None, None))
        if self.verbose:
            if pKa:
                print(f"{self.name} '{description}' pKa={pKa}")
            else:
                print(f"{self.name} does not have a recorded pKa. PKAs = {PKAS}")
        return pKa

    @property
    def pKa_description(self):
        """Return str explaining pKa is relevance (volatile at pH above/below pKa)"""
        pKa, description = PKAS.get(self.name, (None, None))
        return description

    def calc_volatile_portion(self, pH):
        """Return the fraction in the volatile form of mol as a function of pH"""
        pKa, description = PKAS.get(self.name, (None, None))
        HA_to_A_ratio = np.power(10, pKa - pH)
        portion_HA = HA_to_A_ratio / (1 + HA_to_A_ratio)
        if description == "volatile below":
            return portion_HA
        if description == "volatile above":
            return 1 - portion_HA
        raise NotImplementedError(
            f"Don't know how {self.name} can be '{description}' pKa={pKa}"
        )
