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
    def chip_directory(self) -> Path:
        """Get the chip directory."""
        return self.data_directory / "chips"

    @property
    def calibration_directory(self) -> Path:
        """Get the calibration directory."""
        return self.data_directory / "calibrations"

    @property
    def molecule_directory(self) -> Path:
        """Get the molecule director."""
        return self.data_directory / "molecules"

    @property
    def processor_directory(self) -> Path:
        """Get the processor directory."""
        return self.data_directory / "processors"

    @property
    def aux_chip_directory(self) -> Path:
        """Get the chip directory."""
        return self.aux_data_directory / "chips" if self.aux_data_directory else None

    @property
    def aux_calibration_directory(self) -> Path:
        """Get the calibration directory."""
        return self.aux_data_directory / "calibrations" if self.aux_data_directory else None

    @property
    def aux_molecule_directory(self) -> Path:
        """Get the molecule director."""
        return self.aux_data_directory / "molecules" if self.aux_data_directory else None

    @property
    def aux_processor_directory(self) -> Path:
        """Get the processor directory."""
        return self.aux_data_directory / "processors" if self.aux_data_directory else None
