# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains useful shared types for type annotations"""

from pathlib import Path
from typing import TYPE_CHECKING, Dict, Sequence, Union

from .compatability import TypeAlias

if TYPE_CHECKING:
    from spitze.quant import Mixture, Molecule

PATHLIKE: TypeAlias = Union[Path, str]
MOL: TypeAlias = str
MASS: TypeAlias = str
MOLLIST: TypeAlias = Sequence[MOL]
MASSLIST: TypeAlias = Sequence[MASS]
COMPOSITION: TypeAlias = Dict[str, float]
MIXTURE_LIKE: TypeAlias = Union[COMPOSITION, str, "Molecule", "Mixture"]
MASS_TO_SIGNAL: TypeAlias = Dict[MOL, float]
JSONVALUE: TypeAlias = Union[int, float, str, bool]
MOL_TO_FLOAT: TypeAlias = Dict[str, float]
