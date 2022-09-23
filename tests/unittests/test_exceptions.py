# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains unit tests for exceptions.py"""

import inspect
from spectro_inlets_quantification import exceptions


BASE_EXCEPTION = exceptions.QuantError


def test_exceptions_within_hierarchy():
    """Test that all defined exceptions are simple derivatives from the same base"""
    base_arg_spec = inspect.getfullargspec(BASE_EXCEPTION.__init__)
    for name, obj in inspect.getmembers(exceptions):
        if inspect.isclass(obj):
            if obj is BASE_EXCEPTION:
                continue
            assert issubclass(obj, BASE_EXCEPTION)
            # Test that the __init__'s have same signature i.e. not arguments added
            assert inspect.getfullargspec(obj.__init__) == base_arg_spec
