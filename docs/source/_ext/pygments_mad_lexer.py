# -- Custom Lexer ------------------------------------------------------------
from pygments.lexers.scripting import LuaLexer 
from pygments.token import Comment, Operator


class MadLexer(LuaLexer):
    """
    Exactly the same as the Lua Lexer, except for the following changes:
    New name, aliases, filenames and url (also removed mimetypes)
    The comment can also be !
    The character "\\" is now an accepted operator.
    """

    name = 'MAD'
    url = "https://mad.web.cern.ch/mad/"
    aliases = ['mad']
    filenames = ['*.mad']
    mimetypes = []

    tokens = {
        **LuaLexer.tokens,
        'ws': [
            (r'(?:!.*$)', Comment.Single),
            *LuaLexer.tokens["ws"],
        ],
        'base': [
            (r'\\', Operator),
            *LuaLexer.tokens["base"],
        ],
    }