local M = { __author = 'ldeniau', __version = '2015.06.25', help = {}, test = {} }

local ffi = require 'ffi'
local clib = ffi.C
local mabs = math.abs

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

-- Lua inline
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
  if type(y) == 'number' then return complex(x.re / y, x.im / y) end

  local r, d
  if mabs(y.re) < mabs(y.im) then
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
M.__pow = clib.cpow
M.__tostring = M.tostring

complex = ffi.metatype('complex', M)

return complex
