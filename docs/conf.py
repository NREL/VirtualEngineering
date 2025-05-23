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
import os
import sys
sys.path.insert(0, os.path.abspath('../virteng/'))


# -- Project information -----------------------------------------------------

project = 'Virtual Engineering'
copyright = '2021, Ethan Young, Jonathan Stickel, Hari Sitaraman, Andrew Glaws'
author = 'Ethan Young, Jonathan Stickel, Hari Sitaraman, Andrew Glaws'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
			  'sphinx.ext.napoleon',
			  'sphinx.ext.autosectionlabel',
			  'nbsphinx',
			  'jupyter_sphinx',
			  'nbsphinx_link',
			  # 'sphinx_rtd_theme',
			  # 'sphinx.ext.intersphinx',
]

# nbsphinx_execute = 'never'
nbsphinx_requirejs_path = ""
nbsphinx_widgets_path = ""

autodoc_mock_imports = ['numpy', 'scipy', 'fenics-dolfinx', 'virteng', 'ipywidgets', 'pandas']
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.ipynb_checkpoints']
suppress_warnings = ["config.cache"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

html_logo = "figures/ve_logo.png"
html_favicon = "figures/ve_logo_small.png"

autoclass_content = 'both'

