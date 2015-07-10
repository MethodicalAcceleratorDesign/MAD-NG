local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- HELP ------------------------------------------------------------------------

M.__help.self = [[
NAME
  complex

SYNOPSIS
  local complex = require 'complex'
  local I = complex(0,1)
  local z = 2+3i
  local Z = 2+3*I

DESCRIPTION
  The module complex implements the operators and math functions on complex
  numbers:
    (minus) -, +, -, *, /, ^, ==,
    real, imag, conj,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    tostring.

RETURN VALUE
  The constructor of complex numbers

SEE ALSO
  math, gmath
]]

-- DEFS ------------------------------------------------------------------------

local ffi  = require 'ffi'
local clib = require 'mad'

local abs = math.abs

-- Lua API
local complex = ffi.typeof 'complex'
local res = ffi.new 'complex[1]'

function M.real  (x) return x.re end
function M.imag  (x) return x.im end
function M.conj  (x) return complex( x.re, -x.im) end

function M.abs   (x) return clib.mad_cnum_abs (x.re, x.im) end
function M.arg   (x) return clib.mad_cnum_arg (x.re, x.im) end

function M.exp   (x) clib.mad_cnum_exp   (x.re, x.im, res) ; return res[0] end
function M.log   (x) clib.mad_cnum_log   (x.re, x.im, res) ; return res[0] end
function M.sqrt  (x) clib.mad_cnum_sqrt  (x.re, x.im, res) ; return res[0] end
function M.proj  (x) clib.mad_cnum_proj  (x.re, x.im, res) ; return res[0] end

function M.sin   (x) clib.mad_cnum_sin   (x.re, x.im, res) ; return res[0] end
function M.cos   (x) clib.mad_cnum_cos   (x.re, x.im, res) ; return res[0] end
function M.tan   (x) clib.mad_cnum_tan   (x.re, x.im, res) ; return res[0] end
function M.sinh  (x) clib.mad_cnum_sinh  (x.re, x.im, res) ; return res[0] end
function M.cosh  (x) clib.mad_cnum_cosh  (x.re, x.im, res) ; return res[0] end
function M.tanh  (x) clib.mad_cnum_tanh  (x.re, x.im, res) ; return res[0] end

function M.asin  (x) clib.mad_cnum_asin  (x.re, x.im, res) ; return res[0] end
function M.acos  (x) clib.mad_cnum_acos  (x.re, x.im, res) ; return res[0] end
function M.atan  (x) clib.mad_cnum_atan  (x.re, x.im, res) ; return res[0] end
function M.asinh (x) clib.mad_cnum_asinh (x.re, x.im, res) ; return res[0] end
function M.acosh (x) clib.mad_cnum_acosh (x.re, x.im, res) ; return res[0] end
function M.atanh (x) clib.mad_cnum_atanh (x.re, x.im, res) ; return res[0] end

function M.__unm (x)
  return complex(-x.re, -x.im)
end

function M.__eq  (x, y)
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

function M.__pow (x, y)
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

function M.__tostring (x)
      if x.im == 0 then return                        tostring(x.re)
  elseif x.re == 0 then return string.format('%si',                  tostring(x.im))
  elseif x.im <  0 then return string.format('%s%si', tostring(x.re),tostring(x.im))
  else                  return string.format('%s+%si',tostring(x.re),tostring(x.im))
  end
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__index  = M

ffi.metatype('complex', M)

-- END -------------------------------------------------------------------------
return complex
