# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""This module contains unit tests for molecules.py"""
from pathlib import Path
from unittest.mock import patch, mock_open

from attr import fields
from pytest import approx, fixture, mark

from spitze.quant import Medium
from spectro_inlets_quantification.config import Config
from spectro_inlets_quantification.molecule import MoleculeDict, Molecule

MOLECULE_MOD = "spitze.quant.physics.molecule"
CONFIG = Config()
ARGS = (("a", 47), ("b", 42))


@fixture
def molecule():
    """Molecule fixture"""
    return Molecule.load("Ar")


class TestMoleculeDict:
    """Test the MoleculeDict class"""

    def test_init(self):
        """Test the __init__ method"""
        molecule_dict = MoleculeDict(ARGS)
        assert molecule_dict["a"] == 47
        assert set(molecule_dict.keys()) == set(["a", "b"])
        molecule_dict = MoleculeDict(a=47, b=42)
        assert molecule_dict["a"] == 47
        assert isinstance(molecule_dict.medium, Medium)

    def test_get(self):
        """Test get"""
        molecule_dict = MoleculeDict(ARGS)
        assert molecule_dict.get("a") == 47
        with patch.object(Molecule, "load") as mock_load:
            molecule_dict.get("c")
            mock_load.assert_called_once_with("c")

    def test_attrs_post_init(self):
        """Test __attrs_post_init__"""
        molecule = Molecule()
        assert molecule.spectrum is None
        assert molecule.T_of_M is None
        molecule = Molecule(beta=0.8)
        assert molecule.T_of_M is not None
        assert molecule.T_of_M(3.0) == approx(3.0**0.8)
        molecule = Molecule(spectrum_0={"a": 47, "b": 42})
        assert molecule.spectrum == molecule.spectrum_0

    def test_as_dict(self):
        """Test as dict"""
        as_dict = Molecule().as_dict()
        fields_ = fields(Molecule)
        field_names = {
            f.name for f in fields_ if f.init and f.metadata.get("serialize", True)
        }
        field_names == as_dict.keys()

    def test_norm_spectrum(self, molecule):
        """Test calc_norm spectrum"""
        with patch.object(molecule, "calc_norm_spectrum") as mock_calc:
            mock_calc.return_value = 47
            norm_spectrum = molecule.norm_spectrum
            assert norm_spectrum == 47
            mock_calc.assert_called_once()
            norm_spectrum = molecule.norm_spectrum
            assert norm_spectrum == 47
            # Cached property, so still only called once
            mock_calc.assert_called_once()

    @patch("builtins.open", new_callable=mock_open)
    @mark.parametrize("mol_dir", (None, "test_mol_dir"))
    @mark.parametrize("file_name", (None, "test_filename"))
    def test_save(self, mock_open, mol_dir, file_name, molecule):
        """Test save"""
        mol_dir = mol_dir or CONFIG.molecule_directory
        file_name = file_name or molecule.name + ".json"
        save_path = Path(mol_dir) / file_name
        with patch.object(molecule, "as_dict") as mock_as_dict:
            mock_as_dict.return_value = {"a": 47}
            molecule.save(mol_dir=mol_dir, file_name=file_name)
            mock_as_dict.assert_called_once()
            mock_open.assert_called_once_with(save_path, "w")
            handle = mock_open()

            out = "".join(c[0][0] for c in handle.write.call_args_list)
            assert out == '{\n    "a": 47\n}'
