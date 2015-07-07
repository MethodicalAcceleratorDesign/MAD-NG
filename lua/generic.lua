local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  generic -- generic common functions

SYNOPSIS
  local gen = require 'generic'

DESCRIPTION
  The module generic wraps common functions with an object-oriented dispatcher.
  It adds few useful functions:
    isint, inum, round, sign, sinc, step,         (math)
    arg, real, imag, conj, proj,                  (complex) 
    cross, dot                                    (vector)

RETURN VALUES
  The table of generic functions

SEE ALSO
  math
]]

-- DEFS ------------------------------------------------------------------------

local fun = {

-- NUMBER ------------------------------
  -- no wrap
  { 'atan2', 'deg', 'fmod', 'frexp', 'huge', 'isint', 'isnum', 'ldexp', 'max', 'min',
    'modf', 'pi', 'rad', 'random', 'randomseed', 
    fmt = "M.%s = math.%s"
  },
  -- 1 arg
  { 'abs', 'acos', 'asin', 'atan', 'ceil', 'cos', 'cosh', 'exp', 'floor', 'frexp',
    'log10', 'round', 'sign', 'sin', 'sinc', 'sinh', 'sqrt', 'step', 'tan', 'tanh',
    fmt = "local m_%s = math.%s\n" ..
          "M.%s = function (x) return isnum(x) and m_%s(x) or x:%s() end"
  },
  -- 2 args with dispatch on 1st arg (only)
  { 'log', 'pow',
    fmt = "local m_%s = math.%s\n" ..
          "M.%s = function (x, y) return isnum(x) and m_%s(x, y) or x:%s(y) end"
  },

-- OTHER -------------------------------
  { 'arg', 'imag',                          -- complex
    fmt = "M.%s = function (x) return isnum(x) and 0.0 or x:%s() end"
  },

  { 'conj', 'real', 'proj',                 -- complex
    'retain', 'release',                    -- expr
    fmt = "M.%s = function (x) return isnum(x) and x or x:%s() end"
  },

  { 'cross', 'dot',                         -- vector
    fmt = "M.%s = function (x) return x:%s() end"
  },

  { 'tostring',
    fmt = "M.%s = function (x) return isnum(x) and %s(x) or x:%s() end"
  },
}

-- Extensions ------------------------------------------------------------------

local type = type
local abs, ceil, floor, sin = math.abs, math.ceil, math.floor, math.sin
local int_msk = 2^52 + 2^51

local function isnum(x)
  return type(x) == 'number'
end

math.isnum = isnum

math.isint = function (x)
  return isnum(x) and (x + int_msk) - int_msk == x
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

local def = {[[
  local M = ...
  local isnum = math.isnum
]]}

for _, t in ipairs(fun) do
  for _, f in ipairs(t) do
    def[#def+1] = string.gsub(t.fmt, '%%s', f)
  end
end

-- print(table.concat(def, '\n'))
load( table.concat(def, '\n') ) (M)

-- END -------------------------------------------------------------------------
return M
