"""This module contains unit tests for tools.py"""
from time import struct_time
from unittest.mock import patch

import matplotlib.figure
from pytest import mark, approx, raises

from spectro_inlets_quantification.tools import (
    mass_to_M,
    mass_to_pure_mass,
    mass_to_setting,
    tstamp_to_date,
    make_axis,
    Singleton,
)

MASS_DATA = (
    ("M2", 2.0),
    ("M44", 44.0),
    ("M132", 132.0),
    ("M42-", 42.0),
    ("M42-rødgrød", 42.0),
    ("M42A", None),
    ("M32.0", None),
    ("rødgrød", None),
    ("_M32", None),
    ("AM32", None),
)


@mark.parametrize(("value", "result"), MASS_DATA)
def test_mass_to_M(value, result):
    """Test the mass_to_M function"""
    if isinstance(result, float):
        assert mass_to_M(value) == approx(result)
    else:
        with raises(ValueError):
            mass_to_M(value)


@mark.parametrize(("value", "result"), MASS_DATA)
def test_mass_to_pure_mass(value, result):
    """Test the mass_to_pure_mass function"""
    if isinstance(result, float):
        assert mass_to_pure_mass(value) == f"M{int(result)}"
    else:
        with raises(ValueError):
            mass_to_pure_mass(value)


@mark.parametrize(("value", "result"), MASS_DATA)
def test_mass_to_setting(value, result):
    """Test the mass_to_setting function"""
    if isinstance(result, float):
        if "-" in value:
            assert mass_to_setting(value) == value.split("-", maxsplit=1)[1]
        else:
            assert mass_to_setting(value) == "FC"
    else:
        with raises(ValueError):
            mass_to_setting(value)


@mark.parametrize(("month", "month_letter"), ((1, "A"), (7, "G"), (12, "L")))
def test_tstamp_to_date(month, month_letter):
    """Test the tstamp_to_date function"""
    # Mock out time to ensure against local time issues wherever the test runs
    with patch("time.localtime") as mock_localtime:
        mock_localtime.return_value = struct_time(
            (
                2022,
                month,
                4,
                8,
                55,
                36,
                0,
                185,
                1,
            )
        )
        assert tstamp_to_date(1.234567) == f"22{month_letter}04"
    mock_localtime.assert_called_once_with(1.234567)


def test_make_axis():
    """Test the make_axis function"""
    fig, ax = make_axis()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.figure.Axes)


def test_singleton_metaclass():
    """Test the singleton metaclass"""

    class MyClass(metaclass=Singleton):
        @classmethod
        def alternative_constructor(cls):
            return cls()

    assert MyClass() is MyClass()
    assert MyClass() is MyClass.alternative_constructor()
    assert MyClass.alternative_constructor() is MyClass()
