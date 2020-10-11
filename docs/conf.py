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

project = 'DYNAMITE'
copyright = '2020, Stellar Dynamics Group in Vienna'
author = 'Stellar Dynamics Group in Vienna'

# The full version, including alpha/beta/rc tags
release = 'v0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax',
              'sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.inheritance_diagram',
              'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

source_suffix = '.rst'


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
'''
# Uncomment to use alabaster theme
import alabaster
html_theme_path = [alabaster.get_path()]
extensions = ['alabaster']
html_theme = 'alabaster'
'''

html_theme = 'classic'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = 'DynamiteLogo.png'

html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }



# -- Options for LaTeX output --------------------------------------------------

# The paper size ('letter' or 'a4').
#latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
##latex_documents = [
##  ('index', 'sampledoc.tex', u'sampledoc Documentation',
##   u'John Hunter, Fernando Perez, Michael Droettboom', 'manual'),
##]

latex_engine = 'pdflatex'

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
#latex_preamble = ''

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_use_modindex = True

nbsphinx_prolog = """
**Note**:
This page was generated from a Jupyter notebook which can be found at
``docs/{{ env.doc2path(env.docname, base=None) }}``
in your DYNAMITE directory.
----
"""
