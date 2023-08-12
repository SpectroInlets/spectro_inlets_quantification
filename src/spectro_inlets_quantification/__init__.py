# ruff: noqa

# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Import the main classes in spectro_inlets_quantification

Modules in this sub-package have a hierarchy, from lowest to highest level.
A module can not import from a higher-level module.
The order is, from lowest to highest level:

1. constants
2. tools
3. exceptions
4. medium
5. molecule
6. mixture
7. chip
8. peak
9. signal
10. sensitivity
11. calibration
12. quantifier

Several modules and classes are automatically imported with physics, below.

Variables with abbreviated or non-descriptive names, e.g. physical quantities:

* T: Temperature in [K]
* p: pressure in [Pa]
* F(_i)(_M): Mass spec sensitivity in [C/mol]
* F_mat: 2-D matrix of mass spec sensitivity in [C/mol]
* Q_mat: 2-D matrix of inverse mass spec sensitivity in [mol/C]
* H(_i): mass transfer number in [m^3/s]
* c(_i): concentration in [mol/m^3]
* n_dot: flux in [mol/s]
* S(_M): signal in [A]
* sigma(_i): ionization cross section in [Ang^2]                      # not SI!
* f(_i)(_M): predicted sensitivity relative to O2 at M32 [dimensionless]
* M: m/z number in atomic m/z unit [m/z]
* k: constant appearing in equation for f such that f_O2_M32=1 [Ang^-2]
* T_of_M: transmission-amplification function [m/z -> dimensionless]
* alpha: (fitted) predicted sensitivity for F_O2_M32 in [C/mol]
* beta: (fitted) exponent in the transmission function [dimensionless]

Several variables can appear with the subscripts _M and _i

* where _M: mass (e.g. "M32")
* where _i: molecule name (e.g. "O2")

In addition, there are a few abbreviations worth knowing:

* mdict - molecule dictionary, the one and only instance of MoleculeDict
* sf - sensitivity factor, an instance of SensitivityFactor
* sf_list - a list of SensitivityFactor instances
* sensitivity_list - sensitivity list, an instance of SensitivityList
* sm - sensitivity matrix, an instance of SensitivityMatrix
* cal - calibration point, an instance of CalPoint
* cal_list - a list of CalPoint instances

The variable names are from .../Industrial R&D/Quantification/Reports/MS_Theory_v1.0
"""

__version__ = "1.2.3"
__title__ = "spectro-inlets-quantification"
__description__ = "Physics-based tools for quantitative MS"
__author__ = "Spectro Inlets Software"
__email__ = "software@spectroinlets.com"

# I like this to be sure I'm importing from where I think I am:
print(f"importing spectro_inlets_quantification v{__version__} from {__file__}")

# automatic imports, from lowest to highest-level modules:


from . import constants
from . import tools
from . import exceptions

from .medium import Medium
from .molecule import Molecule, MoleculeDict
from .mixture import Mixture, Gas
from .chip import Chip
from .signal import SignalDict, SignalProcessor
from .sensitivity import SensitivityFactor, SensitivityList, SensitivityMatrix
from .calibration import Calibration, CalPoint
from .quantifier import Quantifier
