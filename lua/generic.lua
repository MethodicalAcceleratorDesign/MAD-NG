local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  generic -- generic for common functions

SYNOPSIS
  local gen = require 'generic'

DESCRIPTION
  The module generic wraps common functions with an object-oriented dispatcher.
  It adds few useful functions:
    isint, inum, round, sign, sinc, step,         (math)
    arg, real, imag, conj, proj,                  (complex) 

RETURN VALUES
  The table of generic functions

SEE ALSO
  math, complex
]]

-- DEFS ------------------------------------------------------------------------

local fun = {

-- NUMBERS -----------------------------
  -- no wrap
  { 'unm', 'add', 'sub', 'mul', 'div', 'mod',
    'atan2', 'deg', 'fmod', 'frexp', 'huge', 'isint', 'isnum', 'ldexp', 'max', 'min',
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

-- OBJECTS -----------------------------
  { 'arg', 'imag',                          -- complex
    fmt = "M.%s = function (x) return isnum(x) and 0.0 or x:%s() end"
  },

  { 'conj', 'real', 'proj',                 -- complex
    fmt = "M.%s = function (x) return isnum(x) and x or x:%s() end"
  },

  { 'tostring',
    fmt = "M.%s = function (x) return isnum(x) and %s(x) or x:%s() end"
  },
}

-- Extensions ------------------------------------------------------------------

local int_msk = 2^52 + 2^51

function math.unm   (x)   return -x    end
function math.add   (x,y) return x + y end
function math.sub   (x,y) return x - y end
function math.mul   (x,y) return x * y end
function math.div   (x,y) return x / y end
function math.mod   (x,y) return x % y end

function math.isnum (x)   return type(x) == 'number' end
function math.isint (x)   return math.isnum(x) and (x + int_msk) - int_msk == x end

function math.sign  (x)   return x < 0 and -1 or 1 end
function math.step  (x)   return x < 0 and  0 or 1 end
function math.round (x)   return x < 0 and math.ceil(x-0.5) or math.floor(x+0.5) end
function math.sinc  (x)   return math.abs(x) < 1e-12 and 1.0 or math.sin(x)/x end

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

load( table.concat(def, '\n') ) (M)

-- END -------------------------------------------------------------------------
return M
