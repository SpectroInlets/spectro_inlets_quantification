# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains the configuration object."""

from pathlib import Path
from typing import Optional

from .tools import Singleton

THIS_DIR = Path(__file__).parent


class Config(metaclass=Singleton):
    """Configuration objects."""

    valid_data_file_types = ("calibrations", "chips", "molecules", "processors")

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
    def data_directory(self, path: Optional[Path]) -> None:
        self._data_directory = None if path is None else Path(path)

    @property
    def aux_data_directory(self) -> Path:
        """Get or set data directory."""
        return self._aux_data_directory

    @aux_data_directory.setter
    def aux_data_directory(self, path: Path) -> None:
        self._aux_data_directory = Path(path)

    @property
    def data_directories(self):
        """Return data directories, in the order to look in them."""
        if self.aux_data_directory:
            return [self.aux_data_directory, self.data_directory]
        return [self.data_directory]

    def get_best_data_file(
        self, data_file_type: str, filepath: Path, override_source_dir: Optional[Path] = None
    ) -> Path:
        """Return the best source for a specific data file.

        Args:
            data_file_type (str): One of "calibrations", "chips", "molecules" or "processors"
            filepath (str): The name of the requested data file including extension
            override_source_dir (Path): (Optional) the directory to pull calibrations from
                if given and containing the ``filename``

        'best' means the following order of preferences: [override, aux_dir, packaged] where:

         * **override** is the value of ``override_source_dir``
         * **aux_dir** is self.aux_data_directory / data_file_type
         * **packaged** is self.data_directory / data_file_type

        """
        if data_file_type not in self.valid_data_file_types:
            raise ValueError(
                f"Invalid data file type: '{data_file_type}'. Valid options are: "
                f"{self.valid_data_file_types}"
            )

        sources = [override_source_dir] if override_source_dir else []
        sources += [d / data_file_type for d in self.data_directories]
        for source in sources:
            complete_path = source / filepath
            if complete_path.exists():
                return complete_path

        raise ValueError(
            f"Unable to locate the data file {filepath} of data file type '{data_file_type}' "
            f"in any of the available source directories {sources} for that data file type"
        )

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
