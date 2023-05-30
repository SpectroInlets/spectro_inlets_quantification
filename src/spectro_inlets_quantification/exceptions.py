# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains custom spectro_inlets_quantification exceptions."""


class QuantError(Exception):
    """Base class for exceptions that can arise from SpectroInletsQuant."""


class PeakFitError(QuantError):
    """Error when unable to fit a peak fit (usually due to too low signal)."""


class SensitivityMatrixError(QuantError):
    """Error when unable to invert the sensitivity matrix."""


class MixingError(QuantError):
    """Error when unable to get the capillary flux to match the measurement."""


class RecalibrationError(QuantError):
    """Error when unable to get the sensitivity matrix to match the measurement."""
