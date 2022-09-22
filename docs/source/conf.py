import sys, os, sphinx_rtd_theme
sys.path.append(os.path.abspath("./_ext")) #Add to path here!

# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

master_doc = 'index'
project = 'MAD-NG'
copyright = '2022, Laurent Deniau'
author = 'Laurent Deniau'
release = '0.9.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

primary_domain = "mad" #Use .. default-domain:: c to change to c then .. default-domain:: mad to change back to mad
extensions = ["customRoles", "sphinx-mad-domain"]

source_suffix = {
    '.rst': 'restructuredtext',
    '.mad': 'lua',
}

highlight_language = "lua"

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

html_css_files = [
    'css/custom.css',
]

pygments_style = 'sphinx'

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
  'display_version': True,
  'prev_next_buttons_location': 'both'
}


# -- Options for latexpdf output ----------------------------------------------
latex_toplevel_sectioning = 'part'
latex_elements = {
    'preamble': '\\addto\\captionsenglish{\\renewcommand{\\contentsname}{Table of contents}}',
}

# -- Options for MAN output -------------------------------------------------

man_pages = [
    (master_doc, 'MAD-NG Refence Manual', 'MAD-NG man pages',[author], 1),
    ("sequences", 'Sequence', 'Object man page',[author], 2),
    ("elemfunc", 'Elementary Constants and Functions', 'Elementary Constants and Functions man page',[author], 3),
    #Continually list to get all, could automate this?
]
