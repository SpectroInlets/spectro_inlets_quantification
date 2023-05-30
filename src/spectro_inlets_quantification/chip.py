# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Everything to do with calculating the capillary flux

Variables with abbreviated or non-descriptive names, e.g. physical quantities:

* T: Temperature in [K]
* p: Pressure in [Pa]
* l_cap, w_cap, h_cap: Length, width, and height of capillary in [m]
* n_dot: Capillary flux in [mol/s]
* N_dot: Molecular capillary flux in [s^-1]
* eta: Dynamic viscosity in [Pa*s]
* s: Molecular diameter in [m]
* m: Molecular mass in [kg]

The variable names are from .../Industrial R&D/Quantification/Reports/MS_Theory_v1.0
"""

import json
from math import isclose
from pathlib import Path
from typing import Dict, Union, Optional, Any, Sequence, cast

import numpy as np
from scipy.optimize import fsolve  # type: ignore
from .constants import (
    BOLTZMAN_CONSTANT,
    AVOGADRO_CONSTANT,
    STANDARD_CAPILLARY_WIDTH,
    STANDARD_CAPILLARY_HEIGHT,
    STANDARD_CAPILLARY_LENGTH,
)
from .config import Config
from .exceptions import MixingError
from .medium import Medium
from .mixture import Mixture, Gas
from .custom_types import PATHLIKE, MIXTURE_LIKE, MOL_TO_FLOAT

CONFIG = Config()


class Chip:
    """The Chip class. Mainly used for calculating the capillary flux"""

    # ---- methods whose primary purpose is interface with the .json ---- #
    def __init__(
        self,
        *,
        l_cap: float = STANDARD_CAPILLARY_LENGTH,
        w_cap: float = STANDARD_CAPILLARY_WIDTH,
        h_cap: float = STANDARD_CAPILLARY_HEIGHT,
        carrier: MIXTURE_LIKE = "He",
        solvent: MIXTURE_LIKE = "H2O",
        dry: bool = False,
        verbose: bool = True,
    ) -> None:
        """Create a Chip object given its properties

        Args:
            l_cap (float): Capillary length [m]
            w_cap (float): Capillary width [m]
            h_cap (float): Capillary height [m]
            carrier (str or dict or mixture.Mixture): Carrier gas.
                Strings and dicts are parsed by Mixture.make()
            solvent (str or dict or mixture.Mixture): Solvent.
                Strings and dicts are parsed by Mixture.make()
            dry (bool): Whether the chip is dry (True) or covered in liquid (False)
            verbose (bool): Whether to print stuff to the terminal
        """
        self.verbose = verbose
        self.l_cap = l_cap
        self.w_cap = w_cap
        self.h_cap = h_cap
        self.dry = dry
        self.medium = Medium()
        self._carrier = Gas.make(carrier)
        self._solvent = Mixture.make(solvent)
        self._gas: Optional[Gas] = None

    def __eq__(self, other: object) -> bool:
        """Return whether this Chip is equal to another"""
        if not isinstance(other, self.__class__):
            return False
        for dimension in "lwh":
            dimension_name = f"{dimension}_cap"
            if not isclose(getattr(self, dimension_name), getattr(other, dimension_name)):
                return False
        return True

    def as_dict(self) -> Dict[str, float]:
        """Return a dictionary including everything needed to recreate self"""
        self_as_dict = {"l_cap": self.l_cap}  # more may follow
        return self_as_dict

    def save(self, file_name: str, chip_dir: Optional[PATHLIKE] = None, **kwargs: Any) -> None:
        """Save the `as_dict` form of the chip to a .json file

        Args:
            file_name: Name of the .json file, should include the file suffix. filename.endswith(".json")
            chip_dir: path to directory to save chip in, defaults to :attr:`Config.chip_directory`
                the spitze chips folder.
            kwargs: (other) key word arguments are added to self_as_dict before saving.
        """
        chip_dir = chip_dir or CONFIG.chip_directory
        self_as_dict = self.as_dict()
        self_as_dict.update(kwargs)
        path_to_json = Path(chip_dir) / file_name
        with open(path_to_json, "w") as json_file:
            json.dump(self_as_dict, json_file, indent=4)

    @classmethod
    def load(cls, file_name: str, chip_dir: Optional[PATHLIKE] = None, **kwargs: Any) -> "Chip":
        """loads a chip object from a .json file

        Args:
            file_name: name of the .json file. filename.endswith(".json")
            chip_dir: path to directory to save chip in, defaults to
                the :attr:`Config.chip_directory`
            kwargs: (other) key word arguments are fed to Chip.__init__()

        Returns:
            Chip: a Chip object ready to calculate your capillary flux!
        """
        chip_dir = chip_dir or CONFIG.chip_directory
        path_to_json = (Path(chip_dir) / file_name).with_suffix(".json")
        with open(path_to_json) as json_file:
            self_as_dict = json.load(json_file)
        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    @property
    def carrier(self) -> Gas:
        """dict or str: Carrier gas in the chip. Setting carrier updates mdict."""
        return self._carrier

    @carrier.setter
    def carrier(self, carrier: MIXTURE_LIKE) -> None:
        """Carrier can be set as str or dict, interpreted by `Mixture.make`.

        This sets `chip.gas` to carrier or (if not self.dry) solvent-saturated carrier
        """
        self._carrier = Gas.make(carrier)
        self.reset_gas()

    @property
    def solvent(self) -> Mixture:
        """Mixture: solvent on the chip membrane. Setting solvent updates mdict."""
        return self._solvent

    @solvent.setter
    def solvent(self, solvent: MIXTURE_LIKE) -> None:
        """solvent can be set as str or dict, interpreted by Mixture.make()

        If not self.dry, this sets self.gas to self.carrier saturated with solvent
        """
        self._solvent = Mixture.make(solvent)
        self.reset_gas()

    @property
    def gas(self) -> Gas:
        """Gas: the gas in the chip"""
        if not self._gas:
            self.reset_gas()
        return self._gas

    @gas.setter
    def gas(self, gas: MIXTURE_LIKE) -> None:
        """Gas can be set as str or dict, interpreted by :meth:`Gas.make`"""
        self._gas = Gas.make(gas)

    @property
    def p(self) -> float:
        """Shortcut to ``Chip.medium.p``, as carrier gas is regulated to medium pressure"""
        return self.medium.p

    @p.setter
    def p(self, p: float) -> None:
        """Sets the system pressure in [Pa]."""
        self.medium.p = p

    @property
    def T(self) -> float:
        """Shortcut to ``Chip.medium.T``, since we assume thermal equilibrium"""
        return self.medium.T

    @T.setter
    def T(self, T: float) -> None:
        """Sets the system temperature in [K]."""
        self.medium.T = T

    @property
    def wet(self) -> bool:
        """Whether the chip has liquid on it"""
        return not self.dry

    @wet.setter
    def wet(self, wet: bool) -> None:
        """Set whether the chip has liquid on it"""
        self.dry = not wet

    def reset_gas(self) -> None:
        """Reset chip's gas to carrier or (if not self.dry) solvent-saturated carrier"""
        if self.dry:
            self.gas = self.carrier
        else:
            self.gas = self.carrier.saturated_with(self.solvent, p=self.p, T=self.T)

    # -------------- functions having to do with n_dot ---------------- #

    def calc_n_dot(
        self,
        mol: Optional[str] = None,
        gas: Optional[MIXTURE_LIKE] = None,
        T: Optional[float] = None,
        p: Optional[float] = None,
    ) -> Union[float, MOL_TO_FLOAT]:
        """Calculate the flux through the capillary of a specific molecule

        Args:
            mol (str): The name of a molecule in the chip's gas
            gas (Mixture or dict or str): The gas in the chip, defaults to self.gas
            T (float): Temperature [K], if different from ``medium.T``
            p (float): pressure [Pa], if different from ``medium.p``
        Returns:
            dict: If not mol, return {i: n_dot_i} where n_dot_i is the flux in [mol/s]
            through the capillary of each molecule in the chip's gas
            OR
            float: If mol is given, return the flux in [mol/s] of the specified molecule
        """
        if not gas:
            gas = self.gas
        else:
            gas = Gas.make(gas)
        n_dot_0 = self.calc_n_dot_0(gas=gas, T=T, p=p)
        if mol is None:
            return {mol: x_i * n_dot_0 for mol, x_i in self.gas.comp.items()}
        elif mol in gas.comp:
            return n_dot_0 * gas.comp[mol]
        else:
            print(f"the chip's gas has no {mol}! returning 0.")
            return 0

    def calc_n_dot_0(
        self,
        gas: Optional[MIXTURE_LIKE] = None,
        T: Optional[float] = None,
        p: Optional[float] = None,
    ) -> float:
        """Calculate the total flux through the capillary in [mol/s]

        Uses a weighted average for gas properties and Equation 4.10 of Daniel's Thesis.

        Args:
            gas (dict or str or Mixture): The gas in the chip, defaults to self.gas
            T (float): Temperature [K], if different from ``medium.T``
            p (float): Pressure [Pa], if different from ``medium.p``
        Returns:
            float: The total flux in [mol/s] through the capillary
        """
        N_dot = self.calc_N_dot(T=T, p=p, gas=gas)
        n_dot_0 = N_dot / AVOGADRO_CONSTANT  # converts it from molecules/s to mol/s
        return n_dot_0

    def calc_N_dot(
        self,
        w_cap: Optional[float] = None,
        h_cap: Optional[float] = None,
        l_cap: Optional[float] = None,
        T: Optional[float] = None,
        p: Optional[float] = None,
        gas: Optional[MIXTURE_LIKE] = None,
    ) -> float:
        """Calculate the total molecular flux through the capillary in [s^-1]

        Uses a weighted average for gas properties and Equation 4.10 of Daniel's Thesis.

        Args:
            w_cap (float): Capillary width [m], defaults to self.w_cap
            h_cap (float): Capillary height [m], defaults to self.h_cap
            l_cap (float): Capillary length [m], defaults to self.l_cap
            T (float): Temperature [K], if different from ``medium.T``
            p (float): Pressure [Pa], if different from ``medium.p``
            gas (dict or str): The gas in the chip, defaults to self.gas
        Returns:
            float: The total molecular flux in [s^-1] through the capillary
        """
        if not gas:
            gas_ = self.gas
        else:
            gas_ = Gas.make(gas)

        if w_cap is None:
            w_cap = self.w_cap  # capillary width in [m]
        if h_cap is None:
            h_cap = self.h_cap  # capillary height in [m]
        if l_cap is None:
            l_cap = self.l_cap  # effective capillary length in [m]
        if T is None:
            T = self.T
        if p is None:
            p = self.p
        pi = np.pi
        eta = gas_.eta  # dynamic viscosity in [Pa*s]
        s = gas_.s  # molecule diameter in [m]
        m = gas_.m  # molecule mass in [kg]

        d = ((w_cap * h_cap) / pi) ** 0.5 * 2
        # d = 4.4e-6  #used in Henriksen2009
        a = d / 2
        p_1 = p
        lambda_ = d  # defining the transitional pressure
        # ...from setting mean free path equal to capillary d
        p_t = BOLTZMAN_CONSTANT * T / (2**0.5 * pi * s**2 * lambda_)
        p_2 = 0
        p_m = (p_1 + p_t) / 2  # average pressure in the transitional flow region
        v_m = (8 * BOLTZMAN_CONSTANT * T / (pi * m)) ** 0.5
        # a reciprocal velocity used for short-hand:
        nu = (m / (BOLTZMAN_CONSTANT * T)) ** 0.5

        # ... and now, we're ready for the capillary equation.
        #   (need to turn off black and flake8 for tolerable format)
        # fmt: off
        #   Equation 4.10 of Daniel's Thesis:
        N_dot = (                                                               # noqa
            1 / (BOLTZMAN_CONSTANT * T) * 1 / l_cap * (                         # noqa
                (p_t - p_2) * a**3 * 2 * pi / 3 * v_m + (p_1 - p_t) * (         # noqa
                    a**4 * pi / (8 * eta) * p_m  + a**3 * 2 * pi / 3 * v_m * (  # noqa
                        (1 + 2 * a * nu * p_m / eta) / (                        # noqa
                        1 + 2.48 * a * nu * p_m / eta                           # noqa
                        )                                                       # noqa
                    )                                                           # noqa
                )                                                               # noqa
            )                                                                   # noqa
        )                                                                       # noqa
        # fmt: on
        return cast(float, N_dot)

    # ---- methods for calculating the partial pressures in the chip -------- #

    def calc_pp(
        self,
        n_dot: MOL_TO_FLOAT,
        p: Optional[float] = None,
        T: Optional[float] = None,
        mode: str = "solver",
        relaxed: Sequence[str] = None,
    ) -> MOL_TO_FLOAT:
        """Return the partial pressures of the components of the gas in the chip.

        This method should be called after the flux (n_dot or self.n_dot) has been
        calculated from signals using the sensitivity matrix. It converts the calculated
        fluxes to partial pressures of each molecule in the gas of the chip.
        This can be done in a few ways, specified by "mode":

        naive
            This is the simplest. It assumes that EVERYTHING has a correctly-
            quantified flux. The partial pressure of a molecule in the chip is then just
            the total pressure in the chip times the portion that molecule makes up in
            the sum of the calculated fluxes.

            Advantages:
             * Simple. No use of capillary equation.
             * Robust against uniform (alpha-only) drift of sensitivity.

            Disadvantages:
             * Convolutes uncertainty in each molecule's quantification with the
               uncertainty in the quantification of all other molecules.
             * Fails if there is a significant unquantified or poorly-quantified
               component.

        water
            This normalizes signals according to water's calculated vapour pressure.

            Advantages:
             * Simple. No use of capillary equation.
             * Robust against uniform (alpha-only) drift of sensitivity.

            Disadvantages:
             * Convolutes uncertainty in each molecule's quantification with the
               uncertainty in H2O's quantification, and the calculation of H2O's
               vapour pressure calculation.
             * Fails if there isn't water.

        mix_in
            This is more complex. It lets the flux of the "relaxed" components of
            the carrier gas (typically just He) "float" in order to get the sum of the
            fluxes to match the flux as calculated by the capillary equation.
            See Chip.partial_pressure_by_mix_in for details.

            Advantages:
             * Robust. Works without water, and works okay with significant
               unquantified components (it assumes that they are extra carrier).
             * Uses the sensitivity factors as absolute sensitivity factors - no
               implicit normalization

            Disadvantages:
             * Complex. Uses capillary equation. Requires convergence.
             * Is influenced by uniform sensitivity drift.

        He_solver
            Like mix_in, this adds in the "relaxed" components of the carrier gas
            (typically just He) until the total flux from the MS measurements matches the
            total flux from the capillary calculation. Unlike mix_in, it uses a linear
            solver rather than an iterative approach.

            Advantages:
             * Robust. Works without water, and works okay with significant
               unquantified components (it assumes that they are extra carrier).
             * Uses the sensitivity factors as absolute sensitivity factors - no
               implicit normalization

            Disadvantages:
             * Complex. Uses capillary equation. Depends on solver.
             * Is influenced by uniform sensitivity drift.

        **- - - - - - - - - - End of modes - - - - - - - - - -**

        Args:
            n_dot (dict): {i, n_dot_i}, where n_dot_i is flux in [mol/s] of component i
            p (float): Pressure in [Pa], if different from ``self.medium.p``
            T (float): Temperature in [K], if different from ``self.medium.T``
            mode (str): How to calculate partial pressures.
            relaxed (set of str): Mols to relax if using mix-in mode

        Side Effects:
            ``self.gas`` is reset to match the partial pressures

        Returns
            dict: ``{i: p^i}`` where ``p^i`` is the partial pressure of mol i in [Pa].
        """
        if self.verbose:
            print(f"chip.pp_solver got mode={mode} and n_dot={n_dot}")  # debugging
        if mode == "naive":
            n_dot_total = sum(list(n_dot.values()))
            comp = {mol: n_dot_i / n_dot_total for mol, n_dot_i in n_dot.items()}
            self.gas = Gas.make(comp)
            return self.gas.partial_pressures
        elif mode in ["water", "H2O"]:
            p_H2O = self.gas.mdict["H2O"].calc_p_vap(T=T)
            comp = {mol: p_H2O * n_dot_i / n_dot["H2O"] for mol, n_dot_i in n_dot.items()}
            self.gas = Gas.make(comp)
            return self.gas.partial_pressures
        elif mode in ["mix_in", "mix in"]:
            return self.partial_pressures_by_mix_in(n_dot=n_dot, p=p, T=T, relaxed=relaxed)
        elif mode in ["n_dot_0_solver", "He_solver", "solver"]:
            return self.partial_pressures_by_solver(n_dot=n_dot, p=p, T=T, carrier=self.carrier)
        else:
            raise NotImplementedError(f"calc_pp with mode={mode} is not implemented!")

    def partial_pressures_by_mix_in(
        self,
        n_dot: MOL_TO_FLOAT,
        gas_0: Optional[MIXTURE_LIKE] = None,
        p: Optional[float] = None,
        T: Optional[float] = None,
        relaxed: Sequence[str] = None,
        tolerance: float = 0.01,
        N_loop: float = 20,
        dampening: float = 0.5,
    ) -> MOL_TO_FLOAT:
        """Calculate the composition of the gas in the chip, both analyte and carrier

        This is a difficult, central, and subtle method in the quantification algorithm,
        so it requires some explanation.

        Functions such as ``Calibration.calc_n_dot()`` calculate the flux through the
        chip capillary of analyte molecules based on the MS signals. To calculate
        concentrations, we need to know the mass transfer function, which requires a
        calculation of the total capillary flux. However, the capillary flux depends,
        via the capillary equation, on the composition of the gas in the chip, which is
        some mix of quantified analyte molecules and as-of-yet unquantified (or only
        roughly quantified) molecules in the carrier gas. Naturally added is the requirement
        that the
        total capillary flux is equal to the sum of all the individual fluxes.
        This method takes an iterative approach to solving these almost-circular
        dependencies in calculating the total capillary flux and the composition of the

        The algorithm it uses is:

        0. Start by assuming the gas in the chip only contains carrier gas, and
           define the apparent pressure in the chip to be the measured pressure.
        1. Calculate the total flux (n_dot_0) based on the temperature and
           apparent pressure and the gas assumed to be in the chip. If dampening, mix
           it with the previously calculated n_dot according to dampening.
        2. For each molecule with a flux quantified based on measured MS signals,
           unless that molecule is explicitly "relaxed", assume that molecule's
           fraction in the gas in the chip is equal to its measured flux divided by
           the total flux.
        3. If the sum of the mol fractions is not within the tolerance (default is
           1% of unity):

           a. Set the apparent pressure to the actual pressure times this sum.
           b. Normalize the gas (i.e. divide each mol fraction by this sum), and
           c. Repeat from step (1).

        4. If the sum of the mol fractions of the gas now assumed to be in the chip is
           within the tolerance (default is 1%) of unity, then we're good!

           a. Return the partial pressures

        Steps 1-3 occur in the "mixing loop". The reason that the mixing loop
        converges is as follows: When the analyte fractions are added to the gas
        without normalization (step 2), that will result in a sum of mol fractions
        greater than unity. This will result in an apparent pressure which is higher
        than the actual pressure (step 3), and the total capillary flux calculated in
        the next round (step 1) will be higher. With the same flux of the analyte
        molecules, that results in a smaller fraction for them (step 2),
        and a sum-of-fractions closer to unity, and thus an apparent pressure closer
        to the actual pressure.

        Args:
            n_dot (dict): Trusted fluxes of molecules [pmol/s]. By default, ``n_dot`` will
                be ``self.n_dot_analyte``
            gas_0 (str or dict): Initial guess at the gas in the chip. By default,
                the initial guess is ``self.carrier``
            p (float): Pressure [Pa], if different from ``self.medium.p``
            T (float): Temperature [K], if different from ``self.medium.T``
            relaxed (set of str): Names of the molecules to "relax", see above.
            tolerance (float): How self-consistent the composition of the gas in the
                chip needs to be. Default is 0.01 (i.e. 1%)
            N_loop (int): Maximum number of times in the mixing loop. Default is 20.
            dampening (float): The weighting the previously calculated  ``n_dot_0``
                compared to the new ``n_dot_0``. dampening=0 uses only new. Default is 0.5

        Returns:
            dict: The composition of the gas in the chip with both analyte and carrier.
            On the form ``{i:x^i}`` where ``x^i`` is the mol fraction of molecule i.

        Raises:
            MixingError: If the mixing loop does not converge after ``N_loop`` attempts.
        """

        # parse inputs, evt. getting values from self:
        if p is None:
            p = self.p
        if T is None:
            T = self.T
        n_dot_analyte = n_dot.copy()

        # in this method, we ignore the MS quantification of molecules in self.relax:
        if relaxed:
            for mol in relaxed:
                if mol in n_dot_analyte:
                    n_dot_analyte.pop(mol)

        # Step (0) of the algorithm described in the docstring:
        if gas_0 is None:
            gas = self.gas
        else:
            # copy because we don't want to change the caller's dictionary:
            gas = Gas.make(gas_0)
        # Check which molecules have their flux calculated (only) in this method

        # getting ready for the loop:
        n_loop = 0
        p_app = p
        done = False

        print(f"relaxed = {relaxed}")  # debugging

        # Step (1) of the algorithm described in the doc string
        n_dot_0 = self.calc_n_dot_0(gas=gas, p=p_app, T=T)

        while not done:  # TODO: Optimize convergence / mitigate non-convergence
            print(f"\n\nmixing loop iteration {n_loop}:")  # debugging
            print(f"n_dot_0 = {n_dot_0 * 1e9} nmol/s and \ngas = {gas}")  # debugging

            # Step (2) of the algorithm described in the docstring:
            gas_analyte_comp = {mol: n_dot_i / n_dot_0 for mol, n_dot_i in n_dot_analyte.items()}

            print(f"gas_analyte_comp = {gas_analyte_comp}")  # debugging
            gas.update(gas_analyte_comp)

            total = sum(gas.comp.values())  # should be 1

            if np.isclose(total, 1, atol=tolerance):
                done = True

            # Step (3) of the algorithm described in the docstring:
            gas = Gas.make({mol: val / total for mol, val in gas.comp.items()})
            p_app = p * total
            print(f"apparent pressure = {p_app}, \nactual pressure = {p}")  # debugging

            n_loop += 1
            if n_loop > N_loop:
                raise MixingError("too many times in the mixing loop")

            # re-starting step (1) of the algorithm described in the docstring:
            n_dot_0_new = self.calc_n_dot_0(gas=gas, p=p_app, T=T)
            n_dot_0 = n_dot_0 * dampening + n_dot_0_new * (1 - dampening)

        # Step (4) of the algorithm described in the docstring:
        self.gas = gas
        partial_pressures = {mol: x_i * p for mol, x_i in gas.comp.items()}
        return partial_pressures

    def partial_pressures_by_solver(
        self,
        n_dot: MOL_TO_FLOAT,
        gas_0: Optional[MIXTURE_LIKE] = None,
        p: Optional[float] = None,
        T: Optional[float] = None,
        carrier: Optional[Gas] = None,
    ) -> MOL_TO_FLOAT:
        """Calculate the composition of the gas in the chip, both analyte and carrier

        This is a difficult, central, and subtle method in the quantification algorithm,
        so it requires some explanation.

        Functions such as ``Calibration.calc_n_dot()`` calculate the flux through the
        chip capillary of analyte molecules based on the MS signals. To calculate
        concentrations, we need to know the mass transfer function, which requires a
        calculation of the total capillary flux. However, the capillary flux depends,
        via the capillary equation, on the composition of the gas in the chip, which is
        some mix of quantified analyte molecules and as-of-yet unquantified (or only
        roughly quantified) molecules in the carrier gas. Sanity requires that the
        total capillary flux is equal to the sum of all the individual fluxes.

        Args:
            n_dot (dict): Trusted fluxes of molecules [pmol/s].
            gas_0 (str or dict): Initial guess at the gas in the chip. By default,
                the initial guess is ``self.carrier``
            p (float): Pressure [Pa] if different from ``self.medium.p``
            T (float): Temperature [K] if different from ``self.medium.T``
            carrier (Gas): Gas which to relax

        Returns:
            dict: The composition of the gas in the chip with both analyte and carrier.
            Of the form ``{i:p^i}`` where ``p^i`` is the partial pressure of molecule
            ``i``.

        Raises:
            StopIteration: If the mixing loop does not converge after N_loop attempts.
        """

        # parse inputs, evt. getting values from self:
        if p is None:
            p = self.p
        if T is None:
            T = self.T
        n_dot_analyte = n_dot.copy()

        # in this method, we ignore the (MS-determined) flux of molecules in carrier
        if carrier is None:
            carrier = self.carrier
        for mol in carrier.mol_list:
            if mol in n_dot_analyte:
                n_dot_analyte.pop(mol)
        if gas_0 is None:
            gas_0 = self.gas
        n_dot_total_guess = self.calc_n_dot_0(gas_0, p=p, T=T)

        def make_gas_given_n_dot_total(n_dot_total: float) -> Gas:
            """Return the Gas in the chip given an assumed total capillary flux"""
            comp = {mol: n_dot_i / n_dot_total for mol, n_dot_i in n_dot_analyte.items()}
            mol_fraction_carrier = 1 - sum(comp.values())
            for mol in carrier.mol_list:
                comp[mol] = carrier.comp[mol] * mol_fraction_carrier
            gas = Gas.make(comp)
            return gas

        def n_dot_error(n_dot_total: float) -> float:
            """Return the diff with the implied total flux given an assumed total flux"""
            # fsolve turns the argument into an np.array, causing problems:
            n_dot_total = float(n_dot_total)  # so turn it back into a float
            gas = make_gas_given_n_dot_total(n_dot_total)
            n_dot_0 = self.calc_n_dot_0(gas, T=T, p=p)
            return n_dot_0 - n_dot_total

        n_dot_total = fsolve(n_dot_error, np.array(n_dot_total_guess))[0]
        gas = make_gas_given_n_dot_total(n_dot_total)
        n_dot_0 = self.calc_n_dot_0(gas, p=p, T=T)

        if not np.isclose(n_dot_total, n_dot_0):
            print(f"Warning! n_dot_total={n_dot_total}, but n_dot_0={n_dot_0}")
            raise MixingError

        self.gas = gas
        return {mol: self.p * x_i for mol, x_i in gas.comp.items()}
