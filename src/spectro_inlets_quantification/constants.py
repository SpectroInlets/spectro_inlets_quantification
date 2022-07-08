# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Collection of all relevant constants, to avoid duplicating information"""

# --------- natural constants ------------- #
#: Boltzman constant in [J/K]
BOLTZMAN_CONSTANT = 1.38e-23
#: Avogadro constant in [mol^-1]
AVOGADRO_CONSTANT = 6.0221e23
#: Gas constant in [J/(mol*K)]
GAS_CONSTANT = 8.3143283454
#: This is an empirical dimensionless constant. See :attr:`.Mixture.dynamic_viscosity` and
#: Davidson1993. Davidson writes that 0.375 minimizes error, but goes with 1/3 for simplicity
FLUIDITY_MIXTURE_CONSTANT = 1 / 3

# --------- chip design parameters ------------- #
#: Standard capillary width in [m]
STANDARD_CAPILLARY_WIDTH = 6e-6
#: Standard capillary height in [m]
STANDARD_CAPILLARY_HEIGHT = 6e-6
#: Standard capillary length in [m]
STANDARD_CAPILLARY_LENGTH = 1e-3

# --------- standard external conditions ------------- #
#: Standard temperature in [K]
STANDARD_TEMPERATURE = 298.15
#: Standard temperature in [Pa]
STANDARD_PRESSURE = 1e5
#: Standard vacuum pressure in [Pa]
STANDARD_VACUUM_PRESSURE = 2e-4  # LGACore chamber with 1 bar air through chip [Pa]
#: Standard mixtures
STANDARD_MIXTURES = {
    "air": {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.000412},
    # ^ as of August 8, 2020. See, for changes: https://www.co2.earth/daily-co2
    "cal_gas": {"N2": 0.7, "He": 0.1, "Ar": 0.1, "Xe": 0.1},
    "SCG": {"He": 0.6, "CO2": 0.2, "CH4": 0.2},  # Test Biogas Capillary Gas
}
#: Standard carrier gas
STANDARD_CARRIER_GAS = "He"

# ---------- standard internal conditions ------------- #

STANDARD_SETTING = "FC"  # standard setting specification
"""What this means is FC calibrated so that all peaks are gauss with width=0.7
The m/z=28 peak height in 1 bar air is (20J01) 3.3e-10 [A]
"""
#: Standard ionization energy in [eV]
STANDARD_IONIZATION_ENERGY = 70
#: Standard filament emission current in [A]
STANDARD_IONIZATION_CURRENT = 1e-3
#: Standard CEM voltage in [V]
STANDARD_CEM_VOLTAGE = 1000
#: Stanfardtransmission exponent. See MS_Theory_v1p0.pdf for explanation.
STANDARD_TRANSMISSION_EXPONENT = -1 / 2
#: Standard tuned peak width
STANDARD_WIDTH = 0.7

# ----------- the reference molecule --------------- #
#: The reference molecule
REFERENCE_MOLECULE = "O2"
#: The reference mass
REFERENCE_MASS = "M32"

# ------------- molecule properties ------------------- #
# FIXME: should go in molecule files?

#: pKa values
PKAS = {
    "CO2": (6.3, "volatile below"),
    "NH3": (9.3, "volatile above"),
    "H2S": (7.0, "volatile below"),
    "propionic acid": (4.88, "volatile below"),
    "butyric acid": (4.82, "volatile below"),
    "acetic acid": (4.76, "volatile below"),
}

# ---------- plotting specs ------------------------ #
#: Standard mass plotting colors
STANDARD_COLORS = {
    "M1": "0.8",
    "M2": "b",
    "M4": "m",
    "M18": "y",
    "M28": "0.5",
    "M32": "k",
    "M40": "c",
    "M44": "brown",
    "M15": "r",
    "M26": "g",
    "M27": "limegreen",
    "M30": "darkorange",
    "M31": "yellowgreen",
    "M43": "tan",
    "M45": "darkgreen",
    "M34": "r",
    "M36": "g",
    "M46": "purple",
    "M48": "darkslategray",
    "M20": "slateblue",
    "M16": "steelblue",
    "M19": "teal",
    "M17": "chocolate",
    "M42": "olive",
    "M70": "purple",
    "M3": "orange",
    "M73": "crimson",
    "M74": "r",
    "M60": "g",
    "M58": "darkcyan",
    "M88": "darkred",
    "M89": "darkmagenta",
    "M130": "purple",
    "M132": "purple",
}

#: Standard molecule plotting colors
STANDARD_MOL_COLORS = {
    "O2": "k",
    "N2": "0.6",
    "CO": "0.4",
    "Ar": "c",
    "CO2": "brown",
    "CH4": "r",
    "He": "m",
    "H2": "b",
    "H2O": "y",
    "Xe": "purple",
    "ethanol": "yellowgreen",
    "acetone": "darkkhaki",
    "H2S": "orange",
    "NH3": "purple",
    "propionic acid": "crimson",
    "butyric acid": "teal",
    "acetic acid": "darkslateblue",
}

#: Calibration type plotting markers
CAL_TYPE_SPECS = {
    # calibration type : {kwargs to plt.plot()}. For visualizing all the calibrations
    "initiated": {"marker": "D", "markersize": 10},
    "from mdict": {"marker": "D", "markersize": 10},
    "internal": {"marker": "s", "markersize": 10},
    "semi-internal": {"marker": "^", "markersize": 10},
    "semi": {"marker": "^", "markersize": 10},  # short for semi-internal
    "external": {"marker": "*", "markersize": 12},
    "predicted": {"marker": ".", "markersize": 10},
}
