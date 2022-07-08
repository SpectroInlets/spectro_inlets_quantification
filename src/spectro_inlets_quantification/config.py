# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains the configuration object"""

from pathlib import Path

from .tools import Singleton


THIS_DIR = Path(__file__).parent


class Config(metaclass=Singleton):
    """Configuration objects"""

    def __init__(self, data_directory: Path = THIS_DIR.parent / "data"):
        self._data_directory = Path(data_directory)

    @property
    def data_directory(self) -> Path:
        """Get or set data directory"""
        return self._data_directory

    @data_directory.setter
    def data_directory(self, path: Path) -> None:
        self._data_directory = Path(path)

    @property
    def chip_directory(self) -> Path:
        """Get the chip directory"""
        return self.data_directory / "chips"

    @property
    def calibration_directory(self) -> Path:
        """Get the calibration directory"""
        return self.data_directory / "calibrations"

    @property
    def molecule_directory(self) -> Path:
        """Get the molecule director"""
        return self.data_directory / "molecules"

    @property
    def processor_directory(self) -> Path:
        """Get the processor directory"""
        return self.data_directory / "processors"
