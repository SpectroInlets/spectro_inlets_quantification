"""Boilerplate setup.py

WARNING: This is entirely incomplete and does e.g. not import metadata from __init__ like it
should. It is here only to make the code installable for developers with: pip install -e .

"""
import os
from setuptools import setup, find_packages


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGES = find_packages(where=os.path.join(THIS_DIR, "src"))


setup(
    name="spectro_inlets_quantification",
    version="1.1",
    packages=PACKAGES,
    package_dir={"": "src"},

)
