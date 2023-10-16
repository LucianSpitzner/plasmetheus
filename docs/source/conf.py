# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.append(os.path.abspath('../..'))
# -- Project information

project = 'Plasmetheus'
copyright = '2023, Spitzner'
author = 'Spitzner'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]


autodoc_default_options = {
    'special-members': '__init__'
}

add_module_names = False


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# icon

html_favicon = 'plasIcon.ico'

# -- Options for HTML output

html_theme = 'bootstrap-astropy'

html_theme_options = {
    'logotext1': project,  # white,  semi-bold
    'logotext2': '',  # orange, light
    'logotext3': ':docs'   # white,  light
    }

html_domain_indices = False

# -- Options for EPUB output
epub_show_urls = 'footnote'
