import os
import sys

# ---- Path setup ----
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)

# ---- Project info ----
project = 'Exoplanets_Phase_Curves'
copyright = '2026, Louis-Julien Cartigny'
author = 'Louis-Julien Cartigny'
release = '0.1'

# ---- Extensions ----
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "myst_parser",
]

# ---- Paths ----
templates_path = ['_templates']
exclude_patterns = []

# ---- HTML output ----
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

# ---- Autodoc ----
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}