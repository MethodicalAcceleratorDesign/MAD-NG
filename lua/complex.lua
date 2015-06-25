local M = { __author = 'ldeniau', __version = '2015.06.25', help = {}, test = {} }

local ffi = require 'ffi'
local clib = ffi.C

local cfun = {
  { 'abs', 'arg',
    fmt = "double c%s (double complex x);" },
  {  'exp',  'log',  'proj',  'sqrt',
     'sin',  'cos',  'tan',  'sinh',  'cosh',  'tanh',
    'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh',
    fmt = "double complex c%s (double complex x);" },
  { 'pow',
    fmt = "double complex c%s (double complex x, double complex y);" }
}

-- C API
local cdefs = {}
for _, t in ipairs(cfun) do
  for _, f in ipairs(t) do
    cdefs[#cdefs+1] = t.fmt:format(f)
  end
end
ffi.cdef( table.concat(cdefs, '\n') )

-- Lua API
for _, t in ipairs(cfun) do
  for _, f in ipairs(t) do
    M[f] = clib['c' .. f]
  end
end

local complex

local tocomplex = function (x, y)
  if type(x) == 'number' then x = complex(x) end
  if type(y) == 'number' then y = complex(y) end
  return x, y
end

function M.real (x) return x.re end 
function M.imag (x) return x.im end 
function M.conj (x) return complex(x.re, -x.im) end 

function M.add (x, y)
  x, y = tocomplex(x, y)
  return complex(x.re + y.re, x.im + y.im)
end

function M.sub (x, y)
  x, y = tocomplex(x, y)
  return complex(x.re - y.re, x.im - y.im)
end

function M.mul (x, y)
  x, y = tocomplex(x, y)
  return complex(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re)
end

function M.div (x, y)
  x, y = tocomplex(x, y)
  local r, d
  if math.abs(y.re) < math.abs(y.im) then
    r = y.re / y.im
    d = y.re * r + y.im
    return complex((x.re * r + x.im) / d, (y.im * r - y.re) / d)
  else
    r = y.im / y.re
    d = y.im * r + y.re
    return complex((x.im * r + x.re) / d, (x.im - x.re * r) / d)
  end
end

function M.tostring (x)
-- io.write('complex.__tostring called\n')
      if x.im == 0 then return tostring(x.re)
  elseif x.re == 0 then return string.format('%si',tostring(x.im))
  elseif x.im <  0 then return string.format('%s%si',tostring(x.re),tostring(x.im))
  else                  return string.format('%s+%si',tostring(x.re),tostring(x.im))
  end
end

M.__index = M
M.__eq  = function (x,y) return x.re == y.re and x.im == y.im end
M.__unm = function (x) return complex(-x.re, -x.im) end
M.__add = M.add
M.__sub = M.sub
M.__mul = M.mul
M.__div = M.div
M.__pow = clib.cpow
M.__tostring = M.tostring

complex = ffi.metatype('complex', M)
M.I     = complex(0,1)

return complex
