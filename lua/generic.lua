local M = { __author = 'ldeniau', __version = '2015.06.30', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  generic -- generic common functions

SYNOPSIS
  gen = require 'generic'

DESCRIPTION
  The module generic wraps common functions with a dispatcher

RETURN VALUES
  The table of generic functions

SEE ALSO
  math
]]

-- DEFS ------------------------------------------------------------------------

local fun = {

-- NUMBER ------------------------------
    -- no wrap
    { 'atan2', 'deg', 'fmod', 'frexp', 'huge', 'isint', 'ldexp', 'max', 'min', 'modf', 'pi',
      'rad', 'random', 'randomseed', 
      fmt = [[
          M.%s = math.%s
      ]]
    },
    -- 1 arg
    { 'abs', 'acos', 'asin', 'atan', 'ceil', 'cos', 'cosh', 'exp', 'floor', 'frexp',
      'log10', 'round', 'sign', 'sin', 'sinc', 'sinh', 'sqrt', 'step', 'tan', 'tanh',
      fmt = [[
        local m_%s = math.%s
        M.%s = function (x)
          return type(x) == 'number' and m_%s(x) or x:%s()
        end
      ]]
    },
    -- 2 args with dispatch on 1st arg (only)
    { 'log', 'pow',
      fmt = [[
        local m_%s = math.%s
        M.%s = function (x, y)
          return type(x) == 'number' and m_%s(x, y) or x:%s(y)
        end
      ]]
    },

-- OTHER -------------------------------
    { 'arg', 'imag',                          -- complex
      fmt = [[
        M.%s = function (x) 
          return type(x) == 'number' and 0.0 or x:%s()
        end
      ]]
    },

    { 'conj', 'real', 'proj',                 -- complex
      fmt = [[
        M.%s = function (x) 
          return type(x) == 'number' and x or x:%s()
        end
      ]]
    },

    { 'cross', 'dot',                         -- vector
      fmt = [[
        M.%s = function (x) return x:%s() end
      ]]
    },

    { 'tostring',
      fmt = [[
        M.%s = function (x)
          return type(x) == 'number' and %s(x) or x:%s()
        end
      ]]
    },
}

-- Extensions ------------------------------------------------------------------

local abs, ceil, floor, sin = math.abs, math.ceil, math.floor, math.sin
local int_msk = 2^52 + 2^51

math.isint = function (x)
  return (x + int_msk) - int_msk == x
end

math.round = function (x)
  return x < 0 and ceil(x-0.5) or floor(x+0.5)
end

math.sign = function (x)
  return x < 0 and -1 or 1
end

math.step = function (x)
  return x < 0 and 0 or 1
end

math.sinc = function (x)
  return abs(x) < 1e-12 and 1.0 or sin(x)/x
end

-- Generic API -----------------------------------------------------------------

local def = { "local M = {}" }

for _, t in ipairs(fun) do
  for _, f in ipairs(t) do
    def[#def+1] = string.gsub(t.fmt, '%%s', f)
  end
end

def[#def+1] = "return M"

-- print(table.concat(def, '\n'))

-- END -------------------------------------------------------------------------

local  M = load( table.concat(def, '\n') ) ()
return M
