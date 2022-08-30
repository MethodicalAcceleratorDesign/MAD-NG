from cgitb import text
from unittest import removeResult
from docutils import nodes
from docutils.parsers.rst import directives, Directive
# from sphinx.directives import 

def type_role(name, rawtext, text, lineno, inliner, options=None, context=None):
    node = nodes.emphasis(text = text)
    return [node], []

def macro_role(name, rawtext, text, lineno, inliner, options=None, context=None):
    node = nodes.literal(text = text)
    return [node], []

def unit_role(name, rawtext, text, lineno, inliner, options=None, context=None):
    node = nodes.strong(text = text)
    return [node], []


def setup(app):
    app.add_role("type", type_role)
    app.add_role("macro", macro_role)
    app.add_role("unit", unit_role)

    return
    

