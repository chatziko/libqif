# Sphinx configuration for the qif docs.
#
# The pages are written in markdown (MyST), function signatures come from autodoc (see _ext/qif_docs.py).

import os
import sys

sys.path.insert(0, os.path.abspath('_ext'))

# fail early if the qif module cannot be imported (from PyPI, or in PYTHONPATH when built via cmake), otherwise autodoc silently produces empty docs
import qif

# -- General configuration ------------------------------------------------

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'qif_docs',
]

source_suffix = {'.md': 'markdown', '.rst': 'restructuredtext'}
master_doc = 'index'
exclude_patterns = ['_build', '_notes']
language = 'en'

project = 'qif'
author = 'Kostas Chatzikokolakis'
copyright = author

version = release = qif.__version__     # 'dev' when built from source

# -- MyST -----------------------------------------------------------------

myst_enable_extensions = [
    'dollarmath',       # $...$ and $$...$$ math
    'amsmath',
    'colon_fence',      # ::: directives
    'deflist',
]
myst_heading_anchors = 3

# -- autodoc --------------------------------------------------------------

autodoc_member_order = 'bysource'
add_module_names = False        # posterior(...) instead of qif.measure.bayes_vuln.posterior(...)

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
}

# -- HTML output ----------------------------------------------------------

html_theme = 'furo'
html_title = f'qif {release}'
html_theme_options = {
    'source_repository': 'https://github.com/chatziko/libqif/',
    'source_branch': 'master',
    'source_directory': 'misc/docs/',
}
