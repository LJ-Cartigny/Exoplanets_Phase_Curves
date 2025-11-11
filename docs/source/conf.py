import os
import sys
sys.path.insert(0, os.path.abspath("../../Code_files"))
# sys.path.insert(0, os.path.abspath("../.."))
# print("DEBUG sys.path:", sys.path[:5])  # affiche les premiers chemins pour debug
# print("DEBUG cwd:", os.getcwd())  # affiche le répertoire courant de compilation

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Exoplanets_Phase_Curves'
copyright = '2025, Louis-Julien Cartigny'
author = 'Louis-Julien Cartigny'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    'sphinx.ext.autodoc',
    "sphinx.ext.napoleon", 
    # 'sphinx_autodoc_typehints',
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    # "private-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
    "value": True,
}

autodoc_default_options.update({
    "special-members": "__all__",
    "inherited-members": True,
})


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"
html_theme_options = {"style_nav_header_background": "#2980B9"}
html_static_path = ['_static']

# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_baseurl = "https://lj-cartigny.github.io/Exoplanets_Phase_Curves/"

html_theme_path = []
html_css_files = []

# autodoc_mock_imports = [
#     "numpy",
#     "matplotlib",
#     "pandas",

#     "Code_files.TRAPPIST1_parameters",
# ]

