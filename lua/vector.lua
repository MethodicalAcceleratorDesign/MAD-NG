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
    size, get, set, zeros, ones, set_table,
    real, imag, conj, norm, angle, dot, cross,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    copy, foldl, foldr, map, map2, maps, tostring.

RETURN VALUES
  The constructor of vectors

SEE ALSO
  math, gmath, complex, matrix
]]
 
-- DEFS ------------------------------------------------------------------------

local ffi = require 'ffi'
local gm  = require 'gmath'

-- extend gmath
local istype = ffi.istype
function gm.is_vector  (x) return type(x) == 'cdata' and istype( 'vector_t', x) end
function gm.is_cvector (x) return type(x) == 'cdata' and istype('cvector_t', x) end

-- locals
local isnum, iscpx, iscalar,
      real, imag, conj, ident,
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, add, sub, mul, div, mod, pow = 
      gm.is_number, gm.is_complex, gm.is_scalar,
      gm.real, gm.imag, gm.conj, gm.ident,
      gm.abs, gm.arg, gm.exp, gm.log, gm.sqrt, gm.proj,
      gm.sin, gm.cos, gm.tan, gm.sinh, gm.cosh, gm.tanh,
      gm.asin, gm.acos, gm.atan, gm.asinh, gm.acosh, gm.atanh,
      gm.unm, gm.add, gm.sub, gm.mul, gm.div, gm.mod, gm.pow

local istype, sizeof, fill = ffi.istype, ffi.sizeof, ffi.fill

-- Lua API

ffi.cdef[[
typedef struct { int32_t n; double  data[?]; }  vector_t;
typedef struct { int32_t n; complex data[?]; } cvector_t;
]]

local  vec_ct = ffi.typeof ' vector_t'
local cvec_ct = ffi.typeof 'cvector_t'

local function isvec (x)
  return istype('vector_t', x)
end

local function iscvec (x)
  return istype('cvector_t', x)
end

local function vec_alloc (ct, n)
  local r = ct(n) -- default init: compiled for size <= 128 bytes
  r.n = n
  return r
end

local function vector_from_table(tbl, is_complex)
  local ct = (is_complex or iscpx(tbl[1])) and cvec_ct or vec_ct
  local r = vec_alloc(ct, #tbl)
  return r:set_table(tbl)
end

local function vector (n, is_complex)
  if type(n) == 'table' then
    return vector_from_table(n, is_complex)
  end

  if n > 0 then
    return vec_alloc(is_complex and cvec_ct or vec_ct, n)
  end

  error("invalid argument to vector constructor, expecting size or table")
end

function M.size (x)       return x.n end
function M.get  (x, i)    return x.data[i-1] end
function M.set  (x, i, e) x.data[i-1] = e ; return x end

function M.zeros(x)
  fill(x.data, sizeof(isvec(x) and 'double' or 'complex', x:size()))
  return x
end

function M.ones (x, e)
  e = e or 1 
  for i=0,x:size()-1 do x.data[i] = e end
  return x
end

function M.set_table (x, tbl)
  local n = x:size()
  assert(#tbl == n, "incompatible table size")
  for i=1,n do x.data[i-1] = tbl[i] end
  return x
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
  local n = x:size()
  assert(n == y:size(), "incompatible vector sizes")
  local r = 0
  for i=0,n-1 do r = r + x.data[i] * y.data[i] end
  return r
  -- return M.foldl(M.map2(x, y, gm.mul), 0, gm.add)
end

function M.cross (x, y)
  local n = x:size()
  assert(n == y:size(), "incompatible vector sizes")
  local r = vector(n, iscvec(x) or iscvec(y))
  for i=1,n-2 do
    r.data[i-1] = x.data[i] * y.data[i+1] - x.data[i+1] * y.data[i]
  end
  r.data[n-2] = x.data[n-1] * y.data[0] - x.data[0] * y.data[n-1]
  r.data[n-1] = x.data[  0] * y.data[1] - x.data[1] * y.data[  0]
  return r
end

function M.foldl (x, r, f)
  for i=0,x:size()-1 do r = f(r, x.data[i]) end
  return r
end

function M.foldr (x, r, f)
  for i=0,x:size()-1 do r = f(x.data[i], r) end
  return r
end

function M.map (x, f)
  local n = x:size()
  local r0 = f(x.data[0])
  local r = vector(n, iscpx(r0))
  r.data[0] = r0
  for i=1,n-1 do r.data[i] = f(x.data[i]) end
  return r
end

function M.map2 (x, y, f)
  local n = x:size()
  assert(n == y:size(), "incompatible vector sizes")
  local r0 = f(x.data[0], y.data[0])
  local r = vector(n, iscpx(r0))
  r.data[0] = r0
  for i=1,n-1 do r.data[i] = f(x.data[i], y.data[i]) end
  return r
end

function M.maps (x, y, f)
  if iscalar(x) then
    local n = y:size()    
    local r0 = f(x, y.data[0])
    local r = vector(n, iscpx(r0))
    r.data[0] = r0
    for i=1,n-1 do r.data[i] = f(x, y.data[i]) end
    return r
  end

  if iscalar(y) then
    local n = x:size()    
    local r0 = f(x.data[0], y)
    local r = vector(n, iscpx(r0))
    r.data[0] = r0
    for i=1,n-1 do r.data[i] = f(x.data[i], y) end
    return r
  end

  return M.map2(x, y, f)
end

function M.copy  (x)   return M.map (x,   ident) end
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
  local n = x:size()
  if n ~= y:size() then return false end
  for i=0,n-1 do
    if x.data[i] ~= y.data[i] then return false end
  end
  return true
end

function M.__tostring (x, sep)
  local n = x:size()
  local r = {}
  for i=1,n do r[i] = tostring(x:get(i)) end
  return table.concat(r, sep or ' ')
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__index  = M

ffi.metatype( 'vector_t', M)
ffi.metatype('cvector_t', M)

-- END -------------------------------------------------------------------------
return vector
