import sys, os, sphinx_rtd_theme
sys.path.append(os.path.abspath("./_ext")) #Add to path here!

# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
from pygments_mad_lexer import MadLexer

from sphinx.highlighting import lexers
lexers['mad'] = MadLexer()

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
extensions = [
    "customRoles", 
    "sphinx-mad-domain", 
    ]

source_suffix = {
    '.rst': 'restructuredtext',
}

highlight_language = "mad"
numfig = True

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
latex_docclass = {
    'manual': 'cernrep',
    'howto': 'cernrep'
}
latex_elements = {
    'preamble': '''
    \\usepackage{nameref} % For getting chapter name
    
    \\renewcommand\sphinxtableofcontentshook{}
    \\addto\\captionsenglish{\\renewcommand{\\contentsname}{Table of contents}}

    %table spacing
    %\\setlength{\\tabcolsep}{10pt} % Default value: 6pt
    %\\renewcommand{\\arraystretch}{1.5} % Default value: 1

    %Change the title formats
    \\titleformat{\\chapter      }{\\normalfont\\LARGE\\bfseries}{Chapter \\thechapter . }{1em}{}
    \\titlespacing*{\\chapter}{0pt}{-20pt}{10pt}
    \\titleformat{\\section      }{\\normalfont\\Large\\bfseries}{\\thesection      }{1em}{}
    \\titleformat{\\subsection   }{\\normalfont\\large\\bfseries}{\\thesubsection   }{1em}{}
    \\titleformat{\\subsubsection}{\\normalfont\\large\\bfseries}{\\thesubsubsection}{1em}{}
    
    \\makeatletter
    %Changes headers and footers
    \\fancypagestyle{normal}{ % After page 3 
        \\fancyhf{}
        \\fancyhead[R]{\\thepage}
        \\fancyhead[L]{\\it{\\thechapter . \\MakeUppercase{\\@currentlabelname}}}
        \\renewcommand{\\headrulewidth}{1pt}
        \\renewcommand{\\footrulewidth}{0pt}
    }
    \\fancypagestyle{plain}{ % for up to page 3
        \\fancyhf{}
        \\fancyhead[R]{\\thepage}
    }
    \\makeatother
    ''',
    'fncychap': '',
    'sphinxsetup':  "InnerLinkColor={rgb}{0,0,1}, OuterLinkColor={rgb}{0,0,1}",#'verbatimwithframe = false', #Remove border around code blocks
    'tableofcontents': '\\tableofcontents',
    'maketitle': """
    \\institute{
    Accelerator Beam Physics,\\\\
    CERN, Meyrin, Switzerland.}

    \\begin{abstract}
    The Methodical Accelerator Design -- Next Generation application is an all-in-one standalone versatile tool for particle accelerator design, modeling, and optimization, and for beam dynamics and optics studies. Its general purpose scripting language is based on the simple yet powerful Lua programming language (with a few extensions) and embeds the state-of-art Just-In-Time compiler LuaJIT. Its physics is based on symplectic integration of differential maps made out of GTPSA (Generalized Truncated Power Series). The physics of the transport maps and the normal form analysis were both strongly inspired by the PTC/FPP library from E. Forest. MAD-NG development started in 2016 by the author as a side project of MAD-X, hence MAD-X users should quickly become familiar with its ecosystem, e.g. lattices definition.
    \\begin{center}
    \\texttt{http://cern.ch/mad}
    \\end{center}
    \\end{abstract}

    \\keywords{Methodical Accelerator Design; Accelerator beam physics; Scientific computing; JIT compiler; C and Lua programming.}

    \\maketitle""",
    'extraclassoptions': "oneside"
}
latex_table_style = ['booktabs']
latex_additional_files = ["latex_additional_files/" + x for x in os.listdir("latex_additional_files/")]

# -- Options for MAN output -------------------------------------------------

man_pages = [
    (master_doc, 'MAD-NG Refence Manual', 'MAD-NG man pages',[author], 1),
    ("sequences", 'Sequence', 'Object man page',[author], 2),
    ("elemfunc", 'Elementary Constants and Functions', 'Elementary Constants and Functions man page',[author], 3),
    #Continually list to get all, could automate this?
]
