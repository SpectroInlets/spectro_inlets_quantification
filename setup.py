"""Boilerplate setup.py

WARNING: This is entirely incomplete and does e.g. not import metadata from __init__ like it
should. It is here only to make the code installable for developers with: pip install -e .

"""
import os
from setuptools import setup, find_packages


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGES = find_packages(where=os.path.join(THIS_DIR, "src"))


def read(file):
    """Return the contents of a text file"""
    with open(file) as f:
        return f.read()


setup(
    name="spectro_inlets_quantification",
    # version="1.2",   # this is ignored in favor of that in __init__.py
    packages=PACKAGES,
    package_dir={"": "src"},
    license="DUAL_LICENSE",
    author="Spectro Inlets A/S",
    # author_email=None,
    description="Spectro Inlets open package for MS quantification",
    long_description=read("README.md"),  # note this is overridden by pyproject.toml
    long_description_content_type="text/markdown",
    url="https://github.com/SpectroInlets/spectro_inlets_quantification",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
    ],
    install_requires=read("requirements.txt").split("\n"),
    python_requires=">=3.6",
    include_package_data=True,
    package_data={"": [
        "data/*.yml",
        "data/chips/*.yml",
        "data/molecules/*.yml",
        "data/calibrations/*.yml",
        "data/processors/*.yml"
    ]}
)
