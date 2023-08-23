# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Functional tests for molecule"""
from pathlib import Path

import yaml

from pytest import approx, mark

from spectro_inlets_quantification.config import Config
from spectro_inlets_quantification.molecule import Molecule

CFG = Config()
MOL_NAMES = [
    p.stem
    for p in (CFG.data_directory / "molecules").iterdir()
    if not p.stem.startswith("_") and p.stem != "TEMPLATE"
]
print(MOL_NAMES)


@mark.parametrize("mol_name", MOL_NAMES)
def test_load(mol_name, reset_singletons):
    """Test the load method"""
    file_path = CFG.get_best_data_file("molecules", Path(f"{mol_name}.yml"))
    # Read the raw data
    with open(file_path) as f_:
        raw_data = yaml.safe_load(f_)

    # Sigma needs to have its keys converted back to int, because JSON converts them to strings
    if "sigma" in raw_data:
        raw_data["sigma"] = {int(k): v for k, v in raw_data["sigma"].items()}

    mol = Molecule.load(mol_name)

    for k, v in raw_data.items():
        if isinstance(v, float):
            assert v == approx(getattr(mol, k))
        else:
            assert v == getattr(mol, k)
