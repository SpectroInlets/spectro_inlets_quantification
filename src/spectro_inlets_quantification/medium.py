# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Defines the Medium class.

Medium, unlike the classes of .mixture, is conscious of temperature and pressure

"""
from typing import TYPE_CHECKING, Dict, Optional

from .constants import STANDARD_PRESSURE, STANDARD_TEMPERATURE, STANDARD_VACUUM_PRESSURE
from .tools import Singleton

if TYPE_CHECKING:
    from .custom_types import COMPOSITION
    from .mixture import Mixture
    from .molecule import MoleculeDict

MOLECULE_TO_CONCENTRATION = Dict[str, float]


class Medium(metaclass=Singleton):
    """The liquid on the chip.

    Defines p, T, mixture and concentration for the medium for all of spectro_inlets_quantification

    The class also holds the vacuum side pressure for all of spectro_inlets_quantification.

    Attributes:
        p (float): The medium pressure, defaults to :data:`constants.STANDARD_PRESSURE`
        T (float): The medium temperature, defaults to :data:`.constants.STANDARD_TEMPERATURE`
        p_vac (float): The vacuum pressure, defaults to :data:`constants.STANDARD_VACUUM_PRESSURE`
        mixture (Mixture): The mixture of the medium, defaults to None
        c (dict): The molecule name to concentration mapping, defaults to an empty dict

    """

    def __init__(
        self,
        p: Optional[float] = None,
        T: Optional[float] = None,
        p_vac: Optional[float] = None,
        mixture: Optional["Mixture"] = None,
        c: Optional[MOLECULE_TO_CONCENTRATION] = None,
    ):
        """Initiates the medium. Medium is aware of pressure and temperature.

        It can also be aware of the composition of the medium via mixture and c, though this is
        not yet used.

        Args:
            p (float): Pressure in [Pa]
            T (float): Temperature in [K]
            p_vac (float): The pressure on the vacuum side of the system
            mixture (quant.physics.Mixture): Mixture representing medium's composition
            c (dict): {i: c^i} where c^i is concentration of i in [mM]
        """
        self.p = p or STANDARD_PRESSURE  # the one-and-only system p / [Pa]
        self.T = T or STANDARD_TEMPERATURE  # the one-and-only system T / [K]
        self.p_vac = p_vac or STANDARD_VACUUM_PRESSURE  # the one-and-only system p_vac / [Pa]
        self.mixture = mixture
        if c is None:
            c = {}
        self.c = c

    def __repr__(self) -> str:
        """Return repr string of this object."""
        return (
            f"{self.__class__.__name__}(p={self.p}, T={self.T}, p_vac={self.p_vac}, "
            f"mixture={self.mixture}, c={self.c})"
        )

    @property
    def mdict(self) -> "MoleculeDict":
        """Return the molecule dict."""
        return self.mixture.mdict  # the one-and-only MoleculeDict

    @property
    def comp(self) -> "COMPOSITION":
        """Return the composition."""
        return self.mixture.comp  # the one-and-only relevant composition
