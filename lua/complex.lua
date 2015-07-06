local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

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
     sin,  cos,  tan,  sinh,  cosh,  tanh,
    asin, acos, atan, asinh, acosh, atanh,
    tostring.

RETURN VALUE
  The constructor of complex numbers

SEE ALSO
  math, generic
]]

-- DEFS ------------------------------------------------------------------------

print(package.path)

local ffi  = require 'ffi'
local clib = require 'mad'

local abs = math.abs

-- Lua API
local mt = {}
local res = ffi.new 'complex[1]'
local complex = ffi.typeof 'complex'

function mt.real  (x) return x.re end
function mt.imag  (x) return x.im end
function mt.conj  (x) return complex( x.re, -x.im) end

function mt.abs   (x) return clib.mad_cnum_abs (x.re, x.im) end
function mt.arg   (x) return clib.mad_cnum_arg (x.re, x.im) end

function mt.exp   (x) clib.mad_cnum_exp   (x.re, x.im, res) ; return res[0] end
function mt.log   (x) clib.mad_cnum_log   (x.re, x.im, res) ; return res[0] end
function mt.sqrt  (x) clib.mad_cnum_sqrt  (x.re, x.im, res) ; return res[0] end
function mt.proj  (x) clib.mad_cnum_proj  (x.re, x.im, res) ; return res[0] end

function mt.sin   (x) clib.mad_cnum_sin   (x.re, x.im, res) ; return res[0] end
function mt.cos   (x) clib.mad_cnum_cos   (x.re, x.im, res) ; return res[0] end
function mt.tan   (x) clib.mad_cnum_tan   (x.re, x.im, res) ; return res[0] end
function mt.sinh  (x) clib.mad_cnum_sinh  (x.re, x.im, res) ; return res[0] end
function mt.cosh  (x) clib.mad_cnum_cosh  (x.re, x.im, res) ; return res[0] end
function mt.tanh  (x) clib.mad_cnum_tanh  (x.re, x.im, res) ; return res[0] end

function mt.asin  (x) clib.mad_cnum_asin  (x.re, x.im, res) ; return res[0] end
function mt.acos  (x) clib.mad_cnum_acos  (x.re, x.im, res) ; return res[0] end
function mt.atan  (x) clib.mad_cnum_atan  (x.re, x.im, res) ; return res[0] end
function mt.asinh (x) clib.mad_cnum_asinh (x.re, x.im, res) ; return res[0] end
function mt.acosh (x) clib.mad_cnum_acosh (x.re, x.im, res) ; return res[0] end
function mt.atanh (x) clib.mad_cnum_atanh (x.re, x.im, res) ; return res[0] end

function mt.pow (x, y)
  if type(y) == 'number' then
    if y <  0 then x, y = 1/x, -y end
    if y == 2 then return x*x end
    if y == 1 then return x end
    if y == 0 then return 1 end
    y = complex(y)
  end

  x = complex(x)
  clib.mad_cnum_pow(x.re, x.im, y.re, y.im, res)
  return res[0]
end

function mt.tostring (x)
-- io.write('complex.__tostring called\n') -- __tostring never called (bug?)...
      if x.im == 0 then return                        tostring(x.re)
  elseif x.re == 0 then return string.format('%si',                  tostring(x.im))
  elseif x.im <  0 then return string.format('%s%si', tostring(x.re),tostring(x.im))
  else                  return string.format('%s+%si',tostring(x.re),tostring(x.im))
  end
end

function mt.__unm (x)
  return complex(-x.re, -x.im)
end

function mt.__eq  (x,y)
  x, y = complex(x), complex(y)
  return x.re == y.re and x.im == y.im
end

function mt.__add (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re + y.re, x.im + y.im)
end

function mt.__sub (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re - y.re, x.im - y.im)
end

function mt.__mul (x, y)
  x, y = complex(x), complex(y)
  return complex(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re)
end

function mt.__div (x, y)
  local r, n

  x = complex(x)
  if type(y) == 'number' then
    return complex(x.re / y, x.im / y)
  elseif abs(y.re) < abs(y.im) then
    r = y.re / y.im
    n = 1.0 / (y.re * r + y.im)
    return complex((x.re * r + x.im) * n, (y.im * r - y.re) * n)
  else
    r = y.im / y.re
    n = 1.0 / (y.im * r + y.re)
    return complex((x.im * r + x.re) * n, (x.im - x.re * r) * n)
  end
end

mt.__pow = mt.pow
mt.__tostring = mt.tostring
mt.__index = mt

ffi.metatype('complex', mt)

-- END -------------------------------------------------------------------------
return complex
