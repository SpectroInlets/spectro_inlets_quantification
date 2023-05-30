# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This defines the Mixture and Gas classes

These classes are general tools for dealing with Mixtures

"""

from typing import Dict, Generator, List, Optional, Tuple, Type, TypeVar

import numpy as np

from .constants import (
    FLUIDITY_MIXTURE_CONSTANT,
    STANDARD_MIXTURES,
    STANDARD_PRESSURE,
    STANDARD_TEMPERATURE,
)
from .custom_types import COMPOSITION, MIXTURE_LIKE
from .medium import Medium
from .molecule import Molecule, MoleculeDict

T = TypeVar("T", bound="Mixture")


class Mixture:
    """Class for basic handling of mixtures.

    A powerful wrapper around the composition dictionary comp = {i: x^i}
    """

    def __init__(
        self,
        comp: Optional[COMPOSITION] = None,
        name: Optional[str] = None,
        medium: Optional[Medium] = None,
        verbose: bool = False,
    ) -> None:
        """Initiate a mixture given its composition dictionary and name

        Mixture can get its p and T from self.medium, if medium is given

        Args:
            comp (dict): {i: x^i} where x^i is the mol fraction of molecule i
            name (str): The name of the mixture, defaults to str(comp)
            medium (Medium): the medium from which to read p and T (optional)
            verbose (bool): Whether to print extra info to stdout
        """
        self.comp = comp
        self.name = name or str(comp)
        self.medium = medium
        self.verbose = verbose
        self.mdict = MoleculeDict()  # This is a singleton: all of quant uses same mdict

    @classmethod
    def make(  # noqa: C901
        cls: Type[T],
        mix: MIXTURE_LIKE,
        name: Optional[str] = None,
        verbose: bool = False,
    ) -> T:
        """Return a Mixture object based on mix

        Args:
            mix (dict or str or `Molecule` or `Mixture`): If a dict, it has to be of the
                form {i:x^i} where x^i is the mol fraction of molecule i in the mixture.
                If a str, it should be the name of a standard mixture (a key to
                STANDARD_MIXTURE in the constants module), or the name of a Molecule. If a
                Mixture, it is returned as is.
            name (str): Name to give the mixture, if not derived from mix
            verbose (bool): Whether to print to the terminal when normalizing (debugging use)

        Returns:
            Mixture with the derived comp and name
        """
        if isinstance(mix, cls):
            return mix
        if isinstance(mix, str):
            if name is None:
                name = mix
            if mix in STANDARD_MIXTURES:
                return cls.make(STANDARD_MIXTURES[mix], name=name)
            else:
                comp = {mix: 1}
                return cls(comp=comp, name=name)
        elif isinstance(mix, dict):
            comp = mix.copy()  # don't change the caller's dictionary!
            total = sum(mix.values())
            if not np.isclose(total, 1, rtol=0.01):
                if verbose:
                    print(f"Warning! mix composition added up to {total}. normalizing.")
            for key, val in mix.items():
                comp[key] = val / total
            return cls(comp=comp, name=name)
        else:  # sometimes it's convenient to make a mixture from a Molecule
            try:
                comp = {mix.name: 1}
            except AttributeError:
                if verbose:
                    print(f"Warning!!! Could not make Mixture from {mix}")
                return None
            else:
                if name is None:
                    name = mix.name
                return cls(comp=comp, name=name)

    def __getattr__(self, attr: str) -> float:
        """Unimplemented attributes are weighted average of the Molecule attribute"""
        if "should_not_exist" in attr:
            # Ipython sometimes invokes self._ipython_canary_method_should_not_exist
            # No idea why, but presumably it should raise AttributeError.
            raise AttributeError
        return self.calc_weighted_average(attr)

    def __getitem__(self, key: str) -> Tuple[float, Molecule]:
        """Indexing returns tuple of (fraction, Molecule) for mixture component"""
        return self.comp[key], self.mdict.get(key)

    def __repr__(self) -> str:
        """Return repr string of this object"""
        return f"{self.__class__}({self.name})"

    def __len__(self) -> int:
        """Return length (of self.comp)"""
        return len(self.comp)

    def components(self) -> Generator[Tuple[float, Molecule], None, None]:
        """Iterate over (fraction, molecule) for the components of the mixture"""
        for mol, fraction in self.comp.items():
            yield fraction, self.mdict.get(mol)

    def update(self, comp: COMPOSITION) -> None:
        """Update the mixture with a new composition dictionary"""
        self.comp.update(comp)

    @property
    def mol_list(self) -> List[str]:
        """Return the molecule list"""
        return list(self.comp.keys())

    # T and p can be accessed by the medium, if medium is given
    @property
    def T(self) -> float:
        """Shortcut to Mixture.medium.T. Cannot be set in Mixture"""
        return self.medium.T

    @property
    def p(self) -> float:
        """Shortcut to Mixture.medium.p. Cannot be set in Mixture."""
        return self.medium.p

    # -------- methods which calculate values from the Mixture --------- #

    def calc_weighted_average(self, attr: str) -> float:
        """Estimate a property of a gas as a weighted average of its components

        Molecules without attr are ignored from the average, and if they add up to more
        than 1% of the composition a warning is printed.

        Args:
            attr (str): An attribute of quant.Molecule. The property you want.

        Returns:
            float: A^avg = sum(x^i * A^i), where A is the property and x^i is the mol
            fraction of molecule i in the gas.
        """
        cumulative_attribute = 0.0
        cumulative_fraction = 0.0
        for fraction, molecule in self.components():
            try:
                value = getattr(molecule, attr)
                if value is None:
                    raise AttributeError
                cumulative_attribute += fraction * value
            except AttributeError:
                if self.verbose:
                    print(f"WARNING!!! {molecule.name} has no value for {attr}.")
            else:
                cumulative_fraction += fraction
        if not np.isclose(cumulative_fraction, 1, rtol=0.01) and self.verbose:
            print(f"Warning! Total portion gases with {attr} is {cumulative_fraction}")
        weighted_average = cumulative_attribute / cumulative_fraction
        return weighted_average


class Gas(Mixture):
    """Adds gas-specific methods and attributes to Mixture"""

    def saturated_with(
        self,
        solvent: MIXTURE_LIKE = "H2O",
        p: Optional[float] = None,
        T: Optional[float] = None,
    ) -> "Gas":
        """Return a mixture like self but with the solvent vapor pressure

        If p and T are not given, they come from self.medium. Unless the present gas has
        no medium, in which case p and T are standard conditions.
        This method assumes that solvent is an ideal mixture.

        Args:
            solvent (`Mixture` or mixture-like, see `Mixture.make`): Solvent to saturate with,
                e.g. "H2O", if mixture-like made into Mixture with `Mixture.make`
            p (float): Pressure [Pa]
            T (float): Temperature [K]

        Returns:
            `Gas`: The resulting gas
        """
        if not p:
            p = self.medium.p if self.medium else STANDARD_PRESSURE
        if not T:
            T = self.medium.T if self.medium else STANDARD_TEMPERATURE

        # we aim to get gas_comp, composition of the solvent-saturated gas.
        # it will start as self's comp and than form the solvent's comp.
        solvent = Mixture.make(solvent)
        solvent_comp = solvent.comp
        gas_comp = self.comp.copy()
        p_vap_total = 0  # total vapour pressure / [Pa]

        for mol, fraction in solvent_comp.items():
            molecule = self.mdict.get(mol)
            p_vap = molecule.calc_p_vap(T=T) * fraction  # partial p of mol in [Pa]
            # ^ assumes ideal mixture!
            gas_comp[mol] = p_vap / p  # overwrites if mol was already in gas_comp
            p_vap_total += p_vap
        fraction_dry = 1 - p_vap_total / p  # the portion that is not solvent vapour
        # at this point, gas_comp is not normalized. We need to reduce the amount
        # of each component not in equilibrium with the vapour to normalize it.
        for key, val in gas_comp.items():
            if key not in solvent_comp.keys():
                gas_comp[key] = val * fraction_dry
        # now we're ready to make and return the gas! Let's name it appropriately.
        return Gas(
            gas_comp,
            name=f"{solvent.name}-saturated {self.name}",
            medium=self.medium,
        )

    # -------- methods which calculate values from the Gas --------- #

    @property
    def dynamic_viscosity(self) -> float:
        """The dynamic viscosity in [Pa*s]. Calculated according to Davidson1993

        Davidson, T A. A simple and accurate method for calculating viscosity of gaseous
        mixtures. United States: N. p., 1993. Web.

        And here: https://stacks.cdc.gov/view/cdc/10045/cdc_10045_DS1.pdf

        It uses the fluidity sum equation, Equation 28 of Davidson1993
        It is fully vectorized to improve performance for complex mixtures.
        """
        if len(self.comp) == 1:
            # If only one component, then don't bother with the complex stuff.
            return self.calc_weighted_average("dynamic_viscosity")

        # first, get a list of molecule names and viscosities for molecules with
        # viscosity available:
        mol_list_0 = self.mol_list
        eta_list = []
        mol_list = []
        for mol in mol_list_0:
            eta = self.mdict.get(mol).dynamic_viscosity
            if eta is None:
                if self.verbose:
                    print(
                        f"Removing {mol} from mol_list for viscosity calculation due to "
                        f"viscosity value = {eta}"
                    )
            else:
                mol_list += [mol]
                eta_list += [eta]

        # a vector for the viscosities in [Pa*s] (called mu in Davidson1993):
        eta_vec = np.array(eta_list)
        # a vector for the mol fractions (called x in Davidson1993):
        x_vec = np.array([self.comp[mol] for mol in mol_list])
        # a vector for the masses in [kg] of the mols
        m_vec = np.array([self.mdict.get(mol).m for mol in mol_list])
        # a normalized vector for momentum fraction (Equation 16 of Davidson1993):
        y_vec_not_normalized = x_vec * np.sqrt(m_vec)
        y_vec = y_vec_not_normalized / sum(y_vec_not_normalized)

        # the momentum transfer efficiency matrix (Equation 25 of Davidson1993):
        m_mesh_1, m_mesh_2 = np.meshgrid(m_vec, m_vec)
        efficiency_matrix = 2 * np.sqrt(m_mesh_1 * m_mesh_2) / (m_mesh_1 + m_mesh_2)

        # Equation 28 also needs meshgrids of momentum fraction (y) and viscosity (eta):
        y_mesh_1, y_mesh_2 = np.meshgrid(y_vec, y_vec)
        eta_mesh_1, eta_mesh_2 = np.meshgrid(eta_vec, eta_vec)

        # And here is Equation 28 of Davidson1993 for fluidity:
        fluidity = np.sum(
            np.sum(
                y_mesh_1
                * y_mesh_2
                / (np.sqrt(eta_mesh_1 * eta_mesh_2))
                * np.power(efficiency_matrix, FLUIDITY_MIXTURE_CONSTANT)
            )
        )
        # ... and, viscosity is simply the reciprocal:
        eta = 1 / fluidity
        return eta

    @property
    def eta(self) -> float:
        """Pseudonymn for dynamic_viscosity in [Pa*s]"""
        return self.dynamic_viscosity

    @property
    def partial_pressures(self) -> Dict[str, float]:
        """The partial pressures of the gas's components in [Pa]"""
        return {mol: self.p * x_i for mol, x_i in self.comp.items()}
