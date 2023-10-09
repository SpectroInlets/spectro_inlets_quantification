# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Unit tests for the config module"""
import itertools
from pathlib import Path
from unittest import mock

import pytest

from spectro_inlets_quantification.config import Config


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


def test_get_best_data_file_bad_type():
    config = Config("path")
    with pytest.raises(ValueError):
        config.get_best_data_file("non_existent_type", None)


@pytest.mark.parametrize("exists_return_values", itertools.product((True, False), repeat=3))
def test_get_best_data_file_with_override_and_aux(exists_return_values, reset_singletons):
    # The exists_return_value are all permutations of three booleans, which is used as the
    # return value for whether the generated path exists. Getting all false means the method
    # raises
    data_dir = Path("data_dir")
    aux_dir = Path("aux_dir")
    override_dir = Path("override_dir")
    filename = Path("somefilename.json")

    config = Config(data_dir)
    config.aux_data_directory = aux_dir

    expected = None
    for exists, dirname in zip(
        exists_return_values, (override_dir, aux_dir / "molecules", data_dir / "molecules")
    ):
        if exists:
            expected = dirname / filename
            break

    with mock.patch("pathlib.Path.exists") as exists:
        exists.side_effect = exists_return_values
        if any(exists_return_values):
            return_value = config.get_best_data_file("molecules", filename, override_dir)
            assert return_value == expected
        else:
            with pytest.raises(ValueError):
                config.get_best_data_file("molecules", filename, override_dir)


@pytest.mark.parametrize("exists_return_values", itertools.product((True, False), repeat=2))
def test_get_best_data_file_with_override(exists_return_values, reset_singletons):
    """This test case covers the case where EITHER override or aux is set"""
    # The exists_return_value are all permutations of two booleans, which is used as the
    # return value for whether the generated path exists. Getting all false means the method
    # raises
    data_dir = Path("data_dir")
    override_dir = Path("override_dir")
    filename = Path("somefilename.json")

    config = Config(data_dir)

    expected = None
    for exists, dirname in zip(exists_return_values, (override_dir, data_dir / "molecules")):
        if exists:
            expected = dirname / filename
            break

    with mock.patch("pathlib.Path.exists") as exists:
        exists.side_effect = exists_return_values
        if any(exists_return_values):
            return_value = config.get_best_data_file(
                "molecules", filename, override_source_dir=override_dir
            )
            assert return_value == expected
        else:
            with pytest.raises(ValueError):
                config.get_best_data_file("molecules", filename, override_source_dir=override_dir)


@pytest.mark.parametrize("exists_return_values", itertools.product((True, False), repeat=2))
def test_get_best_data_file_with_aux(exists_return_values, reset_singletons):
    """This test case covers the case where EITHER override or aux is set"""
    # The exists_return_value are all permutations of two booleans, which is used as the
    # return value for whether the generated path exists. Getting all false means the method
    # raises
    data_dir = Path("data_dir")
    aux_dir = Path("aux_dir")
    filename = Path("somefilename.json")

    config = Config(data_dir)

    config.aux_data_directory = aux_dir

    expected = None
    for exists, dirname in zip(
        exists_return_values, (aux_dir / "molecules", data_dir / "molecules")
    ):
        if exists:
            expected = dirname / filename
            break

    with mock.patch("pathlib.Path.exists") as exists:
        exists.side_effect = exists_return_values
        if any(exists_return_values):
            return_value = config.get_best_data_file("molecules", filename)
            assert return_value == expected
        else:
            with pytest.raises(ValueError):
                config.get_best_data_file("molecules", filename)
