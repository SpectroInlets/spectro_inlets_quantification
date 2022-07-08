"""Boilerplate setup.py

WARNING: This is entirely incomplete and does e.g. not import metadata from __init__ like it
should. It is here only to make the code installable for developers with: pip install -e .

"""
import os
from pathlib import Path

from setuptools import setup, find_packages


THIS_DIR = Path(__file__).parent


PACKAGES = find_packages(where=THIS_DIR / "src")

setup(
    name="spectro_inlets_quantification",
    version="1.1",
    packages=PACKAGES,
    package_dir={"": "src"},

)
