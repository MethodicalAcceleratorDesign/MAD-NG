--[=[
 o----------------------------------------------------------------------------o
 |
 | Generic math module
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
  
  Purpose:
  - extends math module and provides object-oriented dispatcher to math
    functions

  Information:
  - defines the 'format' for printing numbers.

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  gmath -- generic math functions

SYNOPSIS
  local gm = require 'gmath'

DESCRIPTION
  The module gmath wraps common math functions with object-oriented dispatch.
  It adds few useful functions:
    is_number, is_integer, is_scalar,                     (math)
    ident, umn, add, sub, mul, div, mod,
    sign, step, round, sinc,
    eq, ne, lt, le, gt, ge,                               (logical)
    arg, real, imag, conj, proj,                          (complex) 
    tostring.                                             (other)

  Other modules may also extend the module gmath (e.g. matrix & tpsa)

RETURN VALUES
  The table of generic functions

SEE ALSO
  math, complex
]]

-- modules -------------------------------------------------------------------o

local ffi = require 'ffi'

-- locals --------------------------------------------------------------------o

local istype = ffi.istype
local format = string.format
local abs, ceil, floor, sin = math.abs, math.ceil, math.floor, math.sin
local int_msk = 2^52 + 2^51

local fun = {

-- Numbers --
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

-- Objects --
  { 'arg', 'imag',                          -- complex
    fmt = "M.%s = function (x) return M.is_number(x) and 0.0 or x:%s() end"
  },

  { 'conj', 'real', 'proj',                 -- complex
    fmt = "M.%s = function (x) return M.is_number(x) and x or x:%s() end"
  },
}

-- implementation ------------------------------------------------------------o

function math.sign  (x) return x < 0 and -1 or 1 end
function math.step  (x) return x < 0 and  0 or 1 end
function math.round (x) return x < 0 and ceil(x-0.5) or floor(x+0.5) end
function math.sinc  (x) return abs(x) < 1e-12 and 1.0 or sin(x)/x end

M.format = "%g" -- default

function M.ident (x)   return  x    end
function M.unm   (x)   return -x    end
function M.add   (x,y) return x + y end
function M.sub   (x,y) return x - y end
function M.mul   (x,y) return x * y end
function M.div   (x,y) return x / y end
function M.mod   (x,y) return x % y end

function M.eq    (x,y) return x == y end
function M.ne    (x,y) return x ~= y end
function M.lt    (x,y) return x <  y end
function M.le    (x,y) return x <= y end
function M.gt    (x,y) return x >  y end
function M.ge    (x,y) return x >= y end

function M.is_number  (x) return type(x) == 'number'   end
function M.is_function(x) return type(x) == 'function' end

function M.is_complex (x) return type(x) == 'cdata' and istype('complex', x) end

function M.is_scalar  (x) return M.is_number(x) or M.is_complex(x) end
function M.is_integer (x) return M.is_number(x) and (x + int_msk) - int_msk == x end

function M.is_table   (x) return type(x) == 'table' and getmetatable(x) == nil end

function M.tostring (x)
  return M.is_number(x) and format(M.format, x) or x.tostring and x:tostring() or tostring(x)
end

-- extend the current module implementation
local def = {[[
  local M = ...
]]}

for _, t in ipairs(fun) do
  for _, f in ipairs(t) do
    def[#def+1] = string.gsub(t.fmt, '%%s', f)
  end
end

load( table.concat(def, '\n') ) (M)

------------------------------------------------------------------------------o
return M
