local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  gmath -- generic math functions

SYNOPSIS
  local gm = require 'gmath'

DESCRIPTION
  The module gmath wraps common math functions with an object-oriented dispatcher.
  It adds few useful functions:
    is_number, is_integer, is_scalar,                     (math)
    ident, umn, add, sub, mul, div, mod,
    sign, step, round, sinc,
    arg, real, imag, conj, proj,                          (complex) 
    tostring.                                             (other)

  Other modules may also extend the module gmath (e.g. vector and matrix)

RETURN VALUES
  The table of generic functions

SEE ALSO
  math, complex, vector, matrix
]]

-- DEFS ------------------------------------------------------------------------

local fun = {

-- NUMBERS -----------------------------
  -- no wrap
  { 'atan2', 'deg', 'fmod', 'frexp', 'huge', 'ldexp', 'max', 'min',
    'modf', 'pi', 'rad', 'random', 'randomseed',
    fmt = "M.%s = math.%s"
  },
  -- 1 arg
  { 'abs', 'acos', 'asin', 'atan', 'ceil', 'cos', 'cosh', 'exp', 'floor', 'frexp',
    'log10', 'round', 'sign', 'sin', 'sinc', 'sinh', 'sqrt', 'step', 'tan', 'tanh',
    fmt = "local m_%s = math.%s\n" ..
          "M.%s = function (x) return M.is_number(x) and m_%s(x) or x:%s() end"
  },
  -- 2 args with dispatch on 1st arg (only)
  { 'log', 'pow',
    fmt = "local m_%s = math.%s\n" ..
          "M.%s = function (x, y) return M.is_number(x) and m_%s(x, y) or x:%s(y) end"
  },

-- OBJECTS -----------------------------
  { 'arg', 'imag',                          -- complex
    fmt = "M.%s = function (x) return M.is_number(x) and 0.0 or x:%s() end"
  },

  { 'conj', 'real', 'proj',                 -- complex
    fmt = "M.%s = function (x) return M.is_number(x) and x or x:%s() end"
  },

  { 'tostring',
    fmt = "M.%s = function (x) return M.is_number(x) and %s(x) or x:%s() end"
  },
}

-- Extensions ------------------------------------------------------------------
local ffi = require 'ffi'

local istype = ffi.istype
local abs, ceil, floor, sin = math.abs, math.ceil, math.floor, math.sin
local int_msk = 2^52 + 2^51

function M.ident (x)   return  x    end
function M.unm   (x)   return -x    end
function M.add   (x,y) return x + y end
function M.sub   (x,y) return x - y end
function M.mul   (x,y) return x * y end
function M.div   (x,y) return x / y end
function M.mod   (x,y) return x % y end

function M.is_number  (x) return type(x) == 'number' end
function M.is_integer (x) return type(x) == 'number' and (x + int_msk) - int_msk == x end
function M.is_complex (x) return type(x) == 'cdata'  and istype('complex', x) end
function M.is_scalar  (x) return M.is_number(x) or M.is_complex(x) end

function math.sign  (x) return x < 0 and -1 or 1 end
function math.step  (x) return x < 0 and  0 or 1 end
function math.round (x) return x < 0 and ceil(x-0.5) or floor(x+0.5) end
function math.sinc  (x) return abs(x) < 1e-12 and 1.0 or sin(x)/x end

-- Generic API -----------------------------------------------------------------

local def = {[[
  local M = ...
]]}

for _, t in ipairs(fun) do
  for _, f in ipairs(t) do
    def[#def+1] = string.gsub(t.fmt, '%%s', f)
  end
end

load( table.concat(def, '\n') ) (M)

-- END -------------------------------------------------------------------------
return M
