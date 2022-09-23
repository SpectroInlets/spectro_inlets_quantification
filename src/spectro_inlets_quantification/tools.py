# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Misc tools

This module contains tools for:
* Matching mass labels
* Making pyplot axis
* Making compact timestamps and the TODAY variable
* And a Singleton metaclass

"""

import re
import time
from typing import Dict, Any, Tuple, TYPE_CHECKING
from math import isclose

if TYPE_CHECKING:
    from matplotlib import pyplot

# Regular expression to match the form "M2", "M44" and "M44-CEM". NOTE: That it is required that
# the first character after the integer mass number, separating the mass from the mass
# spectrometer settings name, be a -

M_MATCH = "^M([0-9]+)([-]{1}.*)?$"


def mass_to_M(mass: str) -> float:
    """Return the m/z value as a float (44.0) given the mass as a str ("M44-CEM")

    Raises:
        ValueError: If mass doesn't adhere to the approved form of either "M44" or "M44-???" where
            ??? is a mass spectrometer settings name e.g. CEM. This name must come after a "-" i.e.
            "M44-CEM"

    """
    match = re.search(M_MATCH, mass)
    if match:
        return float(match.group(1))
    raise ValueError(
        f"Mass '{mass}' does not adhere to approved form e.g: 'M44' or 'M44-???'"
    )


def mass_to_pure_mass(mass: str) -> str:
    """Return just the m/z-specifying mass ("M44") given the full mass ("M44-CEM")

    Raises:
        ValueError: If mass doesn't adhere to the approved form of either "M44" or "M44-???" where
            ??? is a mass spectrometer settings name e.g. CEM. This name must come after a "-" i.e.
            "M44-CEM"

    """
    match = re.match(M_MATCH, mass)
    if match:
        return "M" + match.group(1)
    raise ValueError(
        f"Mass '{mass}' does not adhere to approved form e.g: 'M44' or 'M44-???'"
    )


def mass_to_setting(mass: str) -> str:
    """Return just the mass spectrometer settings name ("CEM") given the full mass ("M44-CEM")"""
    match = re.match(M_MATCH, mass)
    if match:
        if len(match.groups()) > 1 and match.group(2) is not None:
            return match.group(2)[1:]
        else:
            return "FC"
    raise ValueError(
        f"Mass '{mass}' does not adhere to approved form 'M44' or 'M44-???'"
    )


def tstamp_to_date(tstamp: float) -> str:
    """Return the date (str) in special compact form ``yyMdd`` format given the unix time (float)

    In this format the month is given as a capital letter, starting with A for January. E.g. June
    4th, 2022 will become 22G04.

    """
    a = time.localtime(tstamp)
    year = a.tm_year
    month = a.tm_mon
    day = a.tm_mday
    date_string = "{0:02d}{1:1s}{2:02d}".format(
        year % 100, chr(ord("A") + month - 1), day
    )
    return date_string


TODAY = tstamp_to_date(time.time())
"""The date today (i.e. when tools module imported) in compact special form"""


def make_axis() -> Tuple["pyplot.Figure", "pyplot.Axes"]:
    """Make and return a matplotlib.pyplot (figure, axis) pair"""
    from matplotlib import pyplot  # noqa

    return pyplot.subplots()


class Singleton(type):
    """A singleton metaclass"""

    _instances: Dict["Singleton", type] = {}

    def __call__(cls, *args: Any, **kwargs: Any) -> type:
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


def dict_equal_with_close_floats(dict1: Dict[Any, Any], dict2: Dict[Any, Any]) -> bool:
    """Return whether the two dicts are equal, but allow for floats to be only close"""
    if dict1.keys() != dict2.keys():
        return False

    for key, value1 in dict1.items():
        value2 = dict2[key]
        if type(value1) != type(value2):
            return False

        if type(value1) is dict:
            if not dict_equal_with_close_floats(value1, value2):
                return False
        elif type(value1) is list:
            raise RuntimeError(
                "list values not yet supported in `dict_equal_with_close_floats`"
            )
        elif type(value1) is float:
            if not isclose(value1, value2):
                return False
        else:
            if value1 != value2:
                return False

    return True
