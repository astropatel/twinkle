
# list see the documentation:
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

def setup(app):
    app.add_css_file('custom.css')

# nbsphinx_prolog = r"""
# .. raw:: html
#
#     <style>
#         .toctree-wrapper .caption, .toctree-wrapper ul {
#             display: none !important;
#         }
#     </style>
# """
# -- Project information -----------------------------------------------------

project = 'twinkle'
copyright = '2024, Rahul I. Patel'
author = 'Rahul I. Patel'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

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
              ]

templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

html_css_files = ['custom.css']

html_theme_options = {'navigation_depth': 4}

html_logo = '_static/Logo/twinkle_logo_light.png'

# Assuming your `conf.py` has a sibling folder called `_static` with these files
# html_theme_options = {
#    "logo": {
#       "image_light": "_static/Logo/twinkle_logo.png",
#       "image_dark": "_static/Logo/twinkle_logo_invert.png",
#    }
# }