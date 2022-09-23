# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Unit tests for the medium module"""

from unittest.mock import MagicMock, PropertyMock

import pytest

from spectro_inlets_quantification.medium import Medium
from spectro_inlets_quantification.tools import Singleton
from spectro_inlets_quantification.constants import (
    STANDARD_PRESSURE,
    STANDARD_TEMPERATURE,
    STANDARD_VACUUM_PRESSURE,
)

SINGLETONS = (Medium,)


def _reset_singletons():
    """Reset singletons"""
    for singleton in SINGLETONS:
        if singleton in Singleton._instances:
            del Singleton._instances[singleton]


@pytest.fixture(scope="function")
def reset_singletons():
    _reset_singletons()
    yield
    _reset_singletons()


@pytest.mark.parametrize(
    ("name", "default"),
    (
        ("p", STANDARD_PRESSURE),
        ("T", STANDARD_TEMPERATURE),
        ("p_vac", STANDARD_VACUUM_PRESSURE),
        ("mixture", None),
        ("c", {}),
    ),
)
def test_init(reset_singletons, name, default) -> None:
    """Test init"""
    obj = object()
    medium = Medium(**{name: obj})
    assert getattr(medium, name) == obj


@pytest.mark.parametrize(
    ("name", "default"),
    (
        ("p", STANDARD_PRESSURE),
        ("T", STANDARD_TEMPERATURE),
        ("p_vac", STANDARD_VACUUM_PRESSURE),
        ("mixture", None),
        ("c", {}),
    ),
)
def test_init_default_args(reset_singletons, name, default) -> None:
    """Test init"""
    medium = Medium()
    if isinstance(default, float):
        assert getattr(medium, name) == pytest.approx(default)
    else:
        assert getattr(medium, name) == default


def test_mdict_and_comp_properties(reset_singletons) -> None:
    """Test the mdict attribute"""
    mock_mixture = MagicMock()
    mdict_mock = PropertyMock(return_value=47)
    comp_mock = PropertyMock(return_value=42)
    type(mock_mixture).mdict = mdict_mock
    type(mock_mixture).comp = comp_mock
    medium = Medium(mixture=mock_mixture)
    assert medium.mdict == 47
    assert medium.comp == 42
    mdict_mock.assert_called_once_with()
    comp_mock.assert_called_once_with()
