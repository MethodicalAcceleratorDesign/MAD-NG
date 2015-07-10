local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  matrix

SYNOPSIS
  local matrix = require 'matrix'

DESCRIPTION
  The module matrix implements the operators and math functions on matrices:
    (minus) -, +, -, *, /, %, ^, ==, #
    rows, cols, size, sizes, get, set, zeros, ones,
    set_row, set_col, set_diag, set_table,
    real, imag, conj, norm, angle,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    copy, foldl, foldr, map, map2, maps, tostring.

RETURN VALUES
  The constructor of matrices

SEE ALSO
  math, gmath, complex, vector
]]
 
-- DEFS ------------------------------------------------------------------------

local ffi    = require 'ffi'
local vector = require 'vector'
local gm     = require 'gmath'

-- extend gmath
local istype = ffi.istype
function gm.is_matrix  (x) return type(x) == 'cdata' and istype( 'matrix_t', x) end
function gm.is_cmatrix (x) return type(x) == 'cdata' and istype('cmatrix_t', x) end

-- locals
local isnum, iscpx, iscalar, 
      isvec, iscvec, min,
      real, imag, conj, ident,
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, add, sub, mul, div, mod, pow = 
      gm.is_number, gm.is_complex, gm.is_scalar,
      gm.is_vector, gm.is_cvector, gm.min,
      gm.real, gm.imag, gm.conj, gm.ident,
      gm.abs, gm.arg, gm.exp, gm.log, gm.sqrt, gm.proj,
      gm.sin, gm.cos, gm.tan, gm.sinh, gm.cosh, gm.tanh,
      gm.asin, gm.acos, gm.atan, gm.asinh, gm.acosh, gm.atanh,
      gm.unm, gm.add, gm.sub, gm.mul, gm.div, gm.mod, gm.pow

local istype, sizeof, fill = ffi.istype, ffi.sizeof, ffi.fill

-- Lua API

ffi.cdef[[
typedef struct { int32_t nr, nc; double  data[?]; }  matrix_t;
typedef struct { int32_t nr, nc; complex data[?]; } cmatrix_t;
]]

local  mat_ct = ffi.typeof ' matrix_t'
local cmat_ct = ffi.typeof 'cmatrix_t'

local function ismat (x)
  return istype('matrix_t', x)
end

local function iscmat (x)
  return istype('cmatrix_t', x)
end

local function mat_alloc (ct, nr, nc)
  local r = ct(nr*nc) -- default init: compiled for size <= 128 bytes
  r.nr, r.nc = nr, nc
  return r
end

local function matrix_from_table(tbl, is_complex)
  local ct = (is_complex or iscpx(tbl[1][1])) and cmat_ct or mat_ct
  local r = mat_alloc(ct, #tbl, #tbl[1])
  return r:set_table(tbl)
end

local function matrix (nr, nc, is_complex)
  if type(nr) == 'table' then
    return matrix_from_table(nr, is_complex)
  end

  if nr > 0 and nc > 0 then
    return mat_alloc(is_complex and cmat_ct or mat_ct, nr, nc)
  end

  error("invalid argument to matrix constructor, expecting (rows,cols) or table of tables")
end

local function idx (x, i, j)
  return (i-1) * x.nc + (j-1)
end

function M.rows  (x)          return x.nr        end
function M.cols  (x)          return x.nc        end
function M.size  (x)          return x.nr * x.nc end
function M.sizes (x)          return x.nr , x.nc end
function M.get   (x, i, j)    return x.data[idx(x,i,j)] end
function M.set   (x, i, j, e) x.data[idx(x,i,j)] = e ; return x end

function M.zeros(x)
  fill(x.data, sizeof(ismat(x) and 'double' or 'complex', x:size()))
  return x
end

function M.ones (x, e)
  local nr, nc, n = x:sizes()
  x:zeros()
  n, e = min(nr, nc), e or 1
  for i=1,n do x:set(i,i, e) end
  return x
end

function M.set_col (x,jc,v)
  for i=1,x:rows() do x:set(i,jc, v:get(i)) end
  return x
end

function M.set_row (x,ir,v)
  for j=1,x:cols() do x:set(ir,j, v:get(j)) end
  return x
end

function M.set_diag (x,v)
  local nr, nc = x:sizes()
  local n = min(nr, nc)
  for i=1,n do x:set(i,i, v:get(i)) end
  return x
end

function M.set_table (x, tbl)
  local nr, nc = x:sizes()
  for i=1,nr do
    local ti, xi = tbl[i], idx(x,i,0)
    assert(#ti == nc, "incompatible row sizes in table of tables")
    for j=1,nc do x.data[xi + j] = ti[j] end
  end
  return x
end

function M.transpose (x, y)
  local nr, nc = x:sizes()
  local r = y or matrix(nc, nr, iscmat(x))
  for i=1,nr do
    local xi = idx(x,i,0)
    for j=1,nc do r:set(j,i, x.data[xi+j]) end
  end
  return r
end

function M.foldl (x, r, f)
  local nr, nc = x:sizes()
  assert(nr == #r, "incompatible matrix-vector sizes")
  for i=1,nr do
    local ri, xi = r[i], idx(x,i,0)
    for j=1,nc do ri = f(ri, x.data[xi+j]) end
    r[i] = ri
  end
  return r
end

function M.foldr (x, r, f)
  local nr, nc = x:sizes()
  assert(nr == #r, "incompatible matrix-vector sizes")
  for i=1,nr do
    local ri, xi = r[i], idx(x,i,0)
    for j=1,nc do ri = f(x.data[xi+j], ri) end
    r[i] = ri
  end
  return r
end

function M.map (x, f)
  local nr, nc = x:sizes()
  local r0 = f(x.data[0])
  local r = matrix(nr, nc, iscpx(r0))
  r.data[0] = r0
  for i=1,nr*nc-1 do r.data[i] = f(x.data[i]) end
  return r
end

function M.map2 (x, y, f)
  local nr, nc = x:sizes()
  assert(nr == y.nr and nc == y.nc, "incompatible matrix sizes")
  local r0 = f(x.data[0], y.data[0])
  local r = matrix(nr, nc, iscpx(r0))
  r.data[0] = r0
  for i=1,nr*nc-1 do r.data[i] = f(x.data[i], y.data[i]) end
  return r
end

function M.maps (x, y, f)
  if iscalar(x) then
    local nr, nc = y:sizes()    
    local r0 = f(x, y.data[0])
    local r = matrix(nr, nc, iscpx(r0))
    r.data[0] = r0
    for i=1,nr*nc-1 do r.data[i] = f(x, y.data[i]) end
    return r
  end

  if iscalar(y) then
    local nr, nc = x:sizes()    
    local r0 = f(x.data[0], y)
    local r = matrix(nr, nc, iscpx(r0))
    r.data[0] = r0
    for i=1,nr*nc-1 do r.data[i] = f(x.data[i], y) end
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
--function M.__mul (x,y) return M.maps(x,y, mul)   end
--function M.__div (x,y) return M.maps(x,y, div)   end
function M.__mod (x,y) return M.maps(x,y, mod)   end
function M.__pow (x,y) return M.maps(x,y, pow)   end

function M.__eq (x, y)
  local nr, nc = x:sizes()
  if nr ~= y:rows() or nc ~= y:cols() then return false end
  for i=0,nr*nc-1 do
    if x.data[i] ~= y.data[i] then return false end
  end
  return true
end

local function vmul(x,y)
  local n, nr, nc = x:size(), y:sizes()
  assert(n == nr, "incompatible vector-matrix sizes")
  local r = vector(nc, iscvec(x) or iscmat(y))
  for j=1,nc do
    local s, xj = 0, x[j]
    for i=1,nr do s = s + xj * y:get(i,j) end
    r[j] = s
  end
  return r
end

local function mulv(x,y)
  local n, nr, nc = y:size(), x:sizes()
  assert(n == nc, "incompatible matrix-vector sizes")
  local r = vector(nr, iscmat(x) or iscvec(y))
  for i=1,nr do
    local s = 0
    for j=1,nc do s = s + x:get(i,j) * y[j] end
    r[i] = s
  end
  return r
end

local function mmul(x,y)
  local xnr, xnc = x:sizes()
  local ynr, ync = y:sizes()
  assert(xnc == ynr, "incompatible matrices sizes")
  local r = matrix(xnr, ync, iscmat(x) or iscmat(y))
  for i=1,xnr do
    for j=1,ync do
      local s = 0
      for k=1,xnc do s = s + x:get(i,k) * y:get(k,j) end
      r.data[ idx(r,i,j) ] = s
    end
  end
  return r
end

function M.__mul (x,y)
  if iscalar(x) or iscalar(y) then return M.maps(x,y, mul) end
  if isvec(x)   or iscvec(x)  then return   vmul(x,y)      end
  if isvec(y)   or iscvec(y)  then return   mulv(x,y)      end
  return mmul(x,y)
end

function M.__tostring (x, sep, lsep)
  local nr, nc = x:sizes()
  local r, c = {}, {}
  for i=1,nr do
    for j=1,nc do c[j] = tostring(x:get(i,j)) end
    r[i] = table.concat(c, sep or ' ')
  end
  return table.concat(r, lsep or '\n')
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__index  = M

ffi.metatype( 'matrix_t', M)
ffi.metatype('cmatrix_t', M)

-- END -------------------------------------------------------------------------
return matrix
