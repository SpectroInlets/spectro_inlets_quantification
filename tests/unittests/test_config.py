# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Unit tests for the config module"""

from pathlib import Path

import pytest

from spectro_inlets_quantification.config import Config
from spectro_inlets_quantification.tools import Singleton


SINGLETONS = (Config,)


@pytest.fixture(scope="function")
def reset_singletons():
    for singleton in SINGLETONS:
        if singleton in Singleton._instances:
            del Singleton._instances[singleton]
    yield
    for singleton in SINGLETONS:
        if singleton in Singleton._instances:
            del Singleton._instances[singleton]


def test_init(reset_singletons):
    """Test init"""
    path = Path("test")
    config = Config(path)
    assert config._data_directory == path


def test_data_directory_property(reset_singletons):
    """Test data directory property"""
    path = Path("test")
    config = Config(path)
    assert config.data_directory == path
    obj2 = Path("test2")
    config.data_directory = obj2
    assert config.data_directory == obj2


@pytest.mark.parametrize("name", ("chip", "calibration", "molecule", "processor"))
def test_chip_directory(reset_singletons, name) -> None:
    """Test the chip_directory property"""
    path = Path("test")
    config = Config(path)
    assert getattr(config, f"{name}_directory") == path / f"{name}s"
