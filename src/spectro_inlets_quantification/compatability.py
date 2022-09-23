# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains code to ensure compatability over more versions of Python"""

try:
    from typing import TypeAlias
except ImportError:
    from typing_extensions import TypeAlias
