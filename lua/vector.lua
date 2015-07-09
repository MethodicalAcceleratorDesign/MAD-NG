local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  vector

SYNOPSIS
  local vector = require 'vector'

DESCRIPTION
  The module vector implements the operators and math functions on vectors:
    (minus) -, +, -, *, /, %, ^, ==,
    real, imag, conj, norm, angle, dot, cross,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, map, map2, maps, tostring.

RETURN VALUES
  The constructor of vectors

SEE ALSO
  math, complex, generic
]]
 
-- DEFS ------------------------------------------------------------------------

local generic = require 'generic'
local tbl_new = require 'table.new'

local real, imag, conj, 
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, add, sub, mul, div, mod, pow = 
      generic.real, generic.imag, generic.conj,
      generic.abs, generic.arg, generic.exp, generic.log, generic.sqrt, generic.proj,
      generic.sin, generic.cos, generic.tan, generic.sinh, generic.cosh, generic.tanh,
      generic.asin, generic.acos, generic.atan, generic.asinh, generic.acosh, generic.atanh,
      generic.unm, generic.add, generic.sub, generic.mul, generic.div, generic.mod, generic.pow

-- Lua API

local function isvec (x)
  return getmetatable(x) == M
end

local function vector (x)
  if type(x) == 'table'  then
    x[0] = #x
    return setmetatable(x, M)
  end

  if type(x) == 'number' and x > 0 then
    local n = x
    x = tbl_new(n,0)
    x[0] = n
    return setmetatable(x, M)
  end

  error("invalid argument to vector constructor, expecting size or table")
end

function M.norm (x)
  return sqrt(x:dot(x))
end

function M.angle (x, y)
  local w = x:dot(y)
  local v = sqrt(x:dot(x) * y:dot(y))
  return acos(w/v) -- [0, pi]
end

function M.dot (x, y)
  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  local r = 0
  for i=1,n do r = r + x[i] * y[i] end
  return r
  -- return M.foldl(M.map2(x, y, generic.mul), 0, generic.add)
end

function M.cross (x, y)
  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  local r = vector(n)
  for i=1,n-2 do
    r[i] = x[i+1] * y[i+2] - x[i+2] * y[i+1]
  end
  r[n-1] = x[n] * y[1] - x[1] * y[n]
  r[n  ] = x[1] * y[2] - x[2] * y[1]
  return r
end

function M.foldl (x, r, f)
  for i=1,#x do r = f(r, x[i]) end
  return r
end

function M.foldr (x, r, f)
  for i=1,#x do r = f(x[i], r) end
  return r
end

function M.map (x, f)
  local n = #x
  local r = vector(n)
  for i=1,n do r[i] = f(x[i]) end
  return r
end

function M.map2 (x, y, f)
  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  local r = vector(n)
  for i=1,n do r[i] = f(x[i], y[i]) end
  return r
end

function M.maps (x, y, f)
  if not isvec(x) then
    local n = #y
    local r = vector(n)
    for i=1,n do r[i] = f(x, y[i]) end
    return r
  end

  if not isvec(y) then
    local n = #x
    local r = vector(n)
    for i=1,n do r[i] = f(x[i], y) end
    return r
  end

  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  local r = vector(n)
  for i=1,n do r[i] = f(x[i], y[i]) end
  return r
end

function M.real  (x)   return M.map (x,   real)  end
function M.imag  (x)   return M.map (x,   imag)  end
function M.conj  (x)   return M.map (x,   conj)  end
function M.abs   (x)   return M.map (x,   abs)   end
function M.arg   (x)   return M.map (x,   arg)   end
function M.exp   (x)   return M.map (x,   exp)   end
function M.log   (x)   return M.map (x,   log)   end
function M.sqrt  (x)   return M.map (x,   sqrt)  end
function M.proj  (x)   return M.map (x,   proj)  end
function M.sin   (x)   return M.map (x,   sin)   end
function M.cos   (x)   return M.map (x,   cos)   end
function M.tan   (x)   return M.map (x,   tan)   end
function M.sinh  (x)   return M.map (x,   sinh)  end
function M.cosh  (x)   return M.map (x,   cosh)  end
function M.tanh  (x)   return M.map (x,   tanh)  end
function M.asin  (x)   return M.map (x,   asin)  end
function M.acos  (x)   return M.map (x,   acos)  end
function M.atan  (x)   return M.map (x,   atan)  end
function M.asinh (x)   return M.map (x,   asinh) end
function M.acosh (x)   return M.map (x,   acosh) end
function M.atanh (x)   return M.map (x,   atanh) end
function M.__unm (x)   return M.map (x,   unm)   end
function M.__add (x,y) return M.maps(x,y, add)   end
function M.__sub (x,y) return M.maps(x,y, sub)   end
function M.__mul (x,y) return M.maps(x,y, mul)   end
function M.__div (x,y) return M.maps(x,y, div)   end
function M.__mod (x,y) return M.maps(x,y, mod)   end
function M.__pow (x,y) return M.maps(x,y, pow)   end

function M.__eq (x, y)
  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  for i=1,n do
    if x[i] ~= y[i] then return false end
  end
  return true
end

function M.__tostring (x, sep)
  local n = #x
  local r = tbl_new(n,0)
  for i=1,n do r[i] = tostring(x[i]) end
  return table.concat(r, sep or ' ')
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__len    = function (x) return x[0] end
M.__index  = M

-- END -------------------------------------------------------------------------
return vector
