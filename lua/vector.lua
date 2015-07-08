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

local ffi = require 'ffi'
local tbl_new = require 'table.new'
local generic = require 'generic'
local sqrt, acos = generic.sqrt, generic.acos

-- Lua API

local function isvec (x)
  return getmetatable(x) == M
end

local pool = {}

local function vector (x)
  local n = type(x) == 'number' and x or #x
  local ct = pool[n]

  if not ct then
    ffi.cdef("typedef struct { double data[$]; } vector" .. n, n)
    local mt = { __len = n, __index = M }
    ct = ffi.metatype("vector" .. n, mt)
    pool[n] = ct
  end

  local r = type(x) == 'table' and ct({x}) or ct()
end

--[[
local function vector (x)
  if type(x) == 'number' then
    return setmetatable(tbl_new(x,0), M)
  elseif type(x) == 'table' then
    return setmetatable(x, M)
  end
  error("invalid argument to vector constructor, expecting size or table")
end
]]

function M.__eq (x, y)
  local n = #x
  if n ~= #y then error("incompatible vector lengths") end
  for i=1,n do
    if x[i] ~= y[i] then return false end
  end
  return true
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
  local r
  if not isvec(x) then
    local n = #y
    r = vector(n)
    for i=1,n do r[i] = f(x, y[i]) end
  elseif not isvec(y) then
    local n = #x
    r = vector(n)
    for i=1,n do r[i] = f(x[i], y) end
  else
    local n = #x
    r = vector(n)
    if n ~= #y then error("incompatible vector lengths") end
    for i=1,n do r[i] = f(x[i], y[i]) end
  end
  return r
end

function M.real  (x)   return M.map (x,   generic.real)  end
function M.imag  (x)   return M.map (x,   generic.imag)  end
function M.conj  (x)   return M.map (x,   generic.conj)  end

function M.abs   (x)   return M.map (x,   generic.abs)   end
function M.arg   (x)   return M.map (x,   generic.arg)   end
function M.exp   (x)   return M.map (x,   generic.exp)   end
function M.log   (x)   return M.map (x,   generic.log)   end
function M.sqrt  (x)   return M.map (x,   generic.sqrt)  end
function M.proj  (x)   return M.map (x,   generic.proj)  end

function M.sin   (x)   return M.map (x,   generic.sin)   end
function M.cos   (x)   return M.map (x,   generic.cos)   end
function M.tan   (x)   return M.map (x,   generic.tan)   end
function M.sinh  (x)   return M.map (x,   generic.sinh)  end
function M.cosh  (x)   return M.map (x,   generic.cosh)  end
function M.tanh  (x)   return M.map (x,   generic.tanh)  end
function M.asin  (x)   return M.map (x,   generic.asin)  end
function M.acos  (x)   return M.map (x,   generic.acos)  end
function M.atan  (x)   return M.map (x,   generic.atan)  end
function M.asinh (x)   return M.map (x,   generic.asinh) end
function M.acosh (x)   return M.map (x,   generic.acosh) end
function M.atanh (x)   return M.map (x,   generic.atanh) end

function M.__unm (x)   return M.map (x,   generic.unm)   end
function M.__add (x,y) return M.maps(x,y, generic.add)   end
function M.__sub (x,y) return M.maps(x,y, generic.sub)   end
function M.__mul (x,y) return M.maps(x,y, generic.mul)   end
function M.__div (x,y) return M.maps(x,y, generic.div)   end
function M.__mod (x,y) return M.maps(x,y, generic.mod)   end
function M.__pow (x,y) return M.maps(x,y, generic.pow)   end

function M.__tostring (x, sep)
  local n = #x
  local r = tbl_new(n,0)
  for i=1,n do r[i] = tostring(x[i]) end
  return table.concat(r, sep or ' ')
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__index  = M

-- END -------------------------------------------------------------------------
return vector
