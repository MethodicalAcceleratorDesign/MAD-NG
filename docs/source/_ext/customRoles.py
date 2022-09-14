from docutils import nodes
from sphinx import addnodes
from sphinx.locale import _
from sphinx.domains.c import CDomain, CObject, CXRefRole
from sphinx.domains import ObjType


#------May be useful in the future--------
# class CConstObject(CObject):
#     object_type = 'const'

# object_types = CDomain.object_types
# directives = CDomain.directives
# roles = CDomain.roles

# object_types['const'] = ObjType(_('const'),           'const',         'identifier', 'type')
# directives['const'] = CConstObject
# roles['const'] = CXRefRole()
#----------------------------------------

def unit_role(name, rawtext, text, lineno, inliner, options=None, context=None):
    node = nodes.strong(text = text)
    return [node], []

def type_role(name, rawtext, text, lineno, inliner, options=None, context=None):
    node = nodes.emphasis(text = text)
    return [node], []

def setup(app):
    app.add_role("unit", unit_role)
    app.add_role("type", type_role)

    return
    

