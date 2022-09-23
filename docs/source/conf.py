# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "Spectro Inlets Quantification"
copyright = "2022, Spectro Inlets Software Team"
author = "Spectro Inlets Software Team"

# The full version, including alpha/beta/rc tags
release = "1.1"
default_role = "any"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_proof",
]

intersphinx_mapping = {
    # "python": ("https://docs.python.org/3", None),
    # "pandas": ("https://pandas.pydata.org/docs/", None),
    # "numpy": ("https://numpy.org/doc/stable/", None),
    # "matplotlib": ("https://matplotlib.org/stable/", None),
}

# Myst extensions
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    # "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_substitutions = {
    "CH4": "CH$_4$",
    "O2": "O$_2$",
    "C2H4": "C$_2$H$_4$",
    "He": "He",
    "H2": "H$_2$",
    "N2": "N$_2$",
    "Cl2": "Cl$_2$",
    "H2O": "H$_2$O",
    "CO2": "CO_2",
}
suppress_warnings = ["myst.strikethrough"]
numfig = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
