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
# sys.path.append(os.path.abspath('../..'))
# sys.path.append(os.path.abspath('../../twinkle'))
sys.path.insert(0, os.path.abspath('../..'))
# sys.path.insert(0, os.path.abspath('../../twinkle/'))
# sys.path.insert(0, os.path.abspath('../../twinkle/utils'))

# from twinkle import twinkle

def setup(app):
    app.add_css_file('custom.css')


nbsphinx_prolog = r"""
.. raw:: html

    <style>
        .toctree-wrapper .caption, .toctree-wrapper ul {
            display: none !important;
        }
    </style>
"""
# -- Project information -----------------------------------------------------

project = 'twinkle'
copyright = '2024, Rahul I. Patel'
author = 'Rahul I. Patel'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

#mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js'

numfig = True
numfig_secnum_depth = 3

#add napolean extension
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.inheritance_diagram',
              'sphinx.ext.autosummary',
              'sphinx.ext.viewcode',
              'sphinx.builders.html',
              'sphinxcontrib.exceltable',
              'sphinxcontrib.xlsxtable',
              'sphinx.ext.mathjax',
              'nbsphinx']
# 'sphinx.ext.imgmath'
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['custom.css']

html_logo = '_static/Logo/twinkle_logo_light.png'

# Assuming your `conf.py` has a sibling folder called `_static` with these files
html_theme_options = {
   "logo": {
      "image_light": "_static/Logo/twinkle_logo.png",
      "image_dark": "_static/Logo/twinkle_logo_invert.png",
   }
}