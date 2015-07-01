local M = { __author = 'ldeniau', __version = '2015.06.30', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  complex -- complex number

SYNOPSIS
  local complex = require 'complex'
  local I = complex(0,1)
  local z = 2+3i
  local Z = 2+3*I

DESCRIPTION
  The module complex implements the operators and math functions on complex
  numbers:
    (unary) -, (binary) -, +, *, /, ^, ==,
    real, imag, conj,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh,
    tostring.

RETURN VALUES
  The constructor of complex number

SEE ALSO
  math
]]

-- DEFS ------------------------------------------------------------------------

local ffi = require 'ffi'
local bit = require 'bit'
local gen = require 'generic'
local clib = ffi.C
local m_abs, isint, isnum = math.abs, math.isint, math.isnum

local cfun = {
  { 'abs', 'arg',
    fmt = "double c%s (complex);"
  },

  { 'exp',  'log', 'sqrt',  'proj',
    'sin',  'cos',  'tan',  'sinh',  'cosh',  'tanh',
   'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh',
    fmt = "complex c%s (complex);"
  },

  { 'pow',
    fmt = "complex c%s (complex, complex);"
  }
}

-- C API
local cdefs = {}
for _, t in ipairs(cfun) do
  for _, f in ipairs(t) do
    cdefs[#cdefs+1] = t.fmt:format(f)
  end
end
ffi.cdef( table.concat(cdefs, '\n') )

-- C->Lua API
for _, t in ipairs(cfun) do
  for _, f in ipairs(t) do
    M[f] = clib['c' .. f]
  end
end

-- Lua API
local complex

function M.real (x) return x.re end
function M.imag (x) return x.im end
function M.conj (x) return complex(x.re, -x.im) end

function M.__unm (x) return complex(-x.re, -x.im) end

function M.__eq  (x,y)
  x, y = complex(x), complex(y)
  return x.re == y.re and x.im == y.im
end

function M.__add (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re + y.re, x.im + y.im)
end

function M.__sub (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re - y.re, x.im - y.im)
end

function M.__mul (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re)
end

function M.__div (x, y)
  x = complex(x)
  if isnum(y) then
    return complex(x.re / y, x.im / y)
  end

  local r, d
  if m_abs(y.re) < m_abs(y.im) then
    r = y.re / y.im
    d = y.re * r + y.im
    return complex((x.re * r + x.im) / d, (y.im * r - y.re) / d)
  else
    r = y.im / y.re
    d = y.im * r + y.re
    return complex((x.im * r + x.re) / d, (x.im - x.re * r) / d)
  end
end

local band, rshift = bit.band, bit.rshift

function M.pow (x, y)
  if not isint(y) then
    return clib.cpow(x, y)
  end

  if y < 0 then x, y = 1/x, -y end

  local r = 1;

  while true do
    if band(y, 1) ~= 0 then r = r * x end
    y = rshift(y, 1)
    if y == 0 then break end
    x = x * x
  end

  return r
end
M.__pow = M.pow

function M.tostring (x)
-- io.write('complex.__tostring called\n') -- __tostring never called (bug?)...
      if x.im == 0 then return tostring(x.re)
  elseif x.re == 0 then return string.format('%si',tostring(x.im))
  elseif x.im <  0 then return string.format('%s%si',tostring(x.re),tostring(x.im))
  else                  return string.format('%s+%si',tostring(x.re),tostring(x.im))
  end
end

M.__tostring = M.tostring
M.__index = M

-- END -------------------------------------------------------------------------

complex = ffi.metatype('complex', M)

return complex
