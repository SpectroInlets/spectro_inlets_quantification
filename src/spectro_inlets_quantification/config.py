# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains the configuration object."""

from pathlib import Path

from .tools import Singleton

THIS_DIR = Path(__file__).parent


class Config(metaclass=Singleton):
    """Configuration objects."""

    def __init__(self, data_directory: Path = THIS_DIR / "data"):
        """Initialize this objects attributes.

        Args:
            data_directory (Path): The path of the base data directory

        """
        self._data_directory = Path(data_directory)
        self._aux_data_directory = None

    @property
    def data_directory(self) -> Path:
        """Get or set data directory."""
        return self._data_directory

    @data_directory.setter
    def data_directory(self, path: Path) -> None:
        self._data_directory = Path(path)

    @property
    def aux_data_directory(self) -> Path:
        """Get or set data directory."""
        return self._aux_data_directory

    @aux_data_directory.setter
    def aux_data_directory(self, path: Path) -> None:
        self._aux_data_directory = Path(path)

    @property
    def data_directories(self):
        """Return data directories, in the order to look in them"""
        if self.aux_data_directory:
            return [self.aux_data_directory, self.data_directory]
        return [self.data_directory]

    @property
    def chip_directories(self) -> list[Path]:
        """Get the chip directory."""
        return [d / "chips" for d in self.data_directories]

    @property
    def calibration_directories(self) -> list[Path]:
        """Get the calibration directory."""
        return [d / "calibrations" for d in self.data_directories]

    @property
    def molecule_directories(self) -> list[Path]:
        """Get the molecule directory."""
        return [d / "molecules" for d in self.data_directories]

    @property
    def processor_directories(self) -> list[Path]:
        """Get the processor directory."""
        return [d / "processors" for d in self.data_directories]
