local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  matrix

SYNOPSIS
  local matrix = require 'matrix'
  local m1 = matrix {{1,2},{3,4},{5,6}}
  local m2 = matrix(2,3)
  local m3 = m1:transpose()
  m1:transpose(m2) 

DESCRIPTION
  The module matrix implements the operators and math functions on matrices:
    (minus) -, +, -, *, /, %, ^, ==,
    unm, add, sub, mul, div, mod, pow, schur_prod,
    rows, cols, size, sizes, get, set, get0, set0,
    zeros, ones, fill, copy,
    get_row, get_col, get_diag, transpose, t,
    set_row, set_col, set_diag, set_table,
    real, imag, conj, trace, norm, angle,
    dot, inner_prod,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, foreach, map, map2, tostring, totable.

REMARK
  check_bounds  can be set to true to check out of bounds indexes in get , set
  check_bounds0 can be set to true to check out of bounds indexes in get0, set0

RETURN VALUES
  The constructor of matrices

SEE ALSO
  math, gmath, complex, vector, cvector, cmatrix
]]
 
-- DEFS ------------------------------------------------------------------------

local ffi    = require 'ffi'
local linalg = require 'linalg'
local gm     = require 'gmath'

-- locals
local clib            = linalg.cmad
local vector, cvector = linalg.vector, linalg.cvector
local matrix, cmatrix = linalg.matrix, linalg.cmatrix

local isnum, iscpx, iscalar, isvec, iscvec, ismat, iscmat,
      real, imag, conj, ident, min,
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, mod, pow = 
      gm.is_number, gm.is_complex, gm.is_scalar,
      gm.is_vector, gm.is_cvector, gm.is_matrix, gm.is_cmatrix,
      gm.real, gm.imag, gm.conj, gm.ident, gm.min,
      gm.abs, gm.arg, gm.exp, gm.log, gm.sqrt, gm.proj,
      gm.sin, gm.cos, gm.tan, gm.sinh, gm.cosh, gm.tanh,
      gm.asin, gm.acos, gm.atan, gm.asinh, gm.acosh, gm.atanh,
      gm.unm, gm.mod, gm.pow

local istype, cast, sizeof, fill = ffi.istype, ffi.cast, ffi.sizeof, ffi.fill

local cres = ffi.new 'complex[1]'

-- Lua API

local function idx0 (x, i, j)
  if check_bounds0 then
    assert(0 <= i and i < x.nr, "row index out of bounds")
    assert(0 <= j and j < x.nc, "col index out of bounds")
  end
  return i * x.nc + j
end

local function idx1(x, i, j)
  if check_bounds then
    assert(1 <= i and i <= x.nr, "row index out of bounds")
    assert(1 <= j and j <= x.nc, "col index out of bounds")
  end
  return (i-1) * x.nc + (j-1)
end

function M.rows  (x)          return x.nr        end
function M.cols  (x)          return x.nc        end
function M.size  (x)          return x.nr * x.nc end
function M.sizes (x)          return x.nr , x.nc end
function M.get   (x, i, j)    return x.data[idx1(x,i,j)] end
function M.get0  (x, i, j)    return x.data[idx0(x,i,j)] end
function M.set   (x, i, j, e) x.data[idx1(x,i,j)] = e ; return x end
function M.set0  (x, i, j, e) x.data[idx0(x,i,j)] = e ; return x end

function M.reshape (x, nr, nc)
  assert(nr*nc <= x:size(), "incompatible matrix sizes")
  x.nr, x.nc = nr, nc
  return x
end

function M.zeros (x)
  fill(x.data, sizeof(ismat(x) and 'double' or 'complex', x:size()))
  return x
end

function M.ones (x, e_) -- zeros + diagonal
  x:zeros()
  local n, e = min(x:sizes()), e_ or 1
  for i=0,n-1 do x:set0(i,i, e) end
  return x
end

function M.get_col (x, jc, r_)
  local nr = x:rows()
  local r = r_ or ismat(x) and vector(nr) or cvector(nr)
  assert(nr == r:size(), "incompatible matrix-vector sizes")
  for i=1,nr do r:set(i, x:get(i,jc)) end
  return r
end

function M.get_row (x, ir, r_)
  local nc = x:cols()
  local r = r_ or ismat(x) and vector(nr) or cvector(nr)
  assert(nc == r:size(), "incompatible matrix-vector sizes")
  for i=1,nc do r:set(i, x:get(i,ir)) end
  return r
end

function M.get_diag (x, r_)
  local n = min(x:sizes())
  local r = r_ or ismat(x) and vector(n) or cvector(n)
  assert(n == r:size(), "incompatible matrix-vector sizes")
  for i=0,n-1 do r:set0(i, x:get0(i,i)) end
  return r
end

function M.set_col (x, jc, v)
  local nr = x:rows()
  assert(nr == v:size(), "incompatible matrix-vector sizes")
  for i=1,nr do x:set(i,jc, v:get(i)) end
  return x
end

function M.set_row (x, ir, v)
  local nc = x:cols()
  assert(nc == v:size(), "incompatible matrix-vector sizes")
  for j=1,nc do x:set(ir,j, v:get(i)) end
  return x
end

function M.set_diag (x, v)
  local n = min(x:sizes())
  assert(n == v:size(), "incompatible matrix-vector sizes")
  for i=0,n-1 do x:set0(i,i, v:get0(i)) end
  return x
end

function M.set_table (x, tbl)
  local nr, nc = x:sizes()
  assert(#tbl == nr, "incompatible matrix col sizes with table of tables")
  for i=1,nr do
    local ti, xi = tbl[i], idx1(x,i,0)
    assert(#ti == nc, "incompatible matrix row sizes with table of tables")
    for j=1,nc do x.data[xi + j] = ti[j] end
  end
  return x
end

function M.foldl (x, r, f, t_)
  local nr, nc = x:sizes()
  if t_ then
    assert(nc == r:size(), "incompatible matrix-vector sizes")
    for j=0,nc-1 do
      for i=0,nr-1 do r:set0(j, f(r:get0(j), x:get0(i,j))) end
    end
  else
    assert(nr == r:size(), "incompatible matrix-vector sizes")
    for i=0,nr-1 do
      for j=0,nc-1 do r:set0(i, f(r:get0(i), x:get0(i,j))) end
    end
  end
  return r
end

function M.foldr (x, r, f, t_)
  local nr, nc = x:sizes()
  if t_ then
    assert(nc == r:size(), "incompatible matrix-vector sizes")
    for j=0,nc-1 do
      for i=0,nr-1 do r:set0(j, f(x:get0(i,j), r:get0(j))) end
    end
  else
    assert(nr == r:size(), "incompatible matrix-vector sizes")
    for i=0,nr-1 do
      for j=0,nc-1 do r:set0(i, f(x:get0(i,j), r:get0(i))) end
    end
  end
  return r
end

function M.foreach (x, f)
  local nr, nc = x:sizes()
  for i=1,nr do
    for i=1,nc do f(x:get(i,j), i, j) end
  end
  return x
end

function M.map (x, f, r_)
  local nr, nc = x:sizes()
  local r0 = f(x:get(1,1))
  local r = r_ or isnum(r0) and matrix(nr, nc) or cmatrix(nr, nc)
  assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
  r.data[0] = r0
  for i=1,nr*nc-1 do r.data[i] = f(x.data[i]) end
  return r
end

function M.map2 (x, y, f, r_)
  assert(ismat(y) or iscmat(y), "matrix expected for 2nd argument")
  local nr, nc = x:sizes()
  local r0 = f(x:get(1,1), y:get(1,1))
  local r = r_ or isnum(r0) and matrix(nr, nc) or cmatrix(nr, nc)
  assert(nr == y:rows() and nc == y:cols(), "incompatible matrix sizes")
  assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
  r.data[0] = r0
  for i=1,nr*nc-1 do r.data[i] = f(x.data[i], y.data[i]) end
  return r
end

function M.transpose (x, r_)
  local nr, nc = x:sizes()
  local r = r_ or ismat(x) and matrix(nc, nr) or cmatrix(nc, nr)
  assert(nr == r:cols() and nc == r:rows(), "incompatible matrix sizes")
  if ismat(x) then
    clib.mad_mat_trans (x.data, r.data, nr, nc)
  else
    clib.mad_cmat_trans(x.data, r.data, nr, nc)
  end
  return r
end

M.t = M.transpose

function M.trace (x)
  local n = min(x:sizes())
  local r = 0
  for i=0,n-1 do r = r + x:get0(i,i) end
  return r
end

function M.inner_prod (x, y)
  return (x:transpose() * y):trace()
end

M.dot = M.inner_prod

function M.norm (x)
  -- Frobenius (consistent with inner_prod)
  if ismat(x) then
    return sqrt(clib.mad_vec_dot(x.data, x.data, x:size()))
  else
    clib.mad_cvec_dot(x.data, x.data, cres, x:size())
    return sqrt(cres[0])
  end
end

function M.angle (x, y)
  local w = x:inner_prod(y)
  local v = x:norm() * y:norm()
  return acos(w / v) -- [0, pi]
end

function M.fill  (x, e )  return x:map(function () return e end, x) end
function M.copy  (x, r_)  return x:map(ident, r_) end
function M.real  (x, r_)  return x:map(real , r_) end
function M.imag  (x, r_)  return x:map(imag , r_) end
function M.conj  (x, r_)  return x:map(conj , r_) end
function M.abs   (x, r_)  return x:map(abs  , r_) end
function M.arg   (x, r_)  return x:map(arg  , r_) end
function M.exp   (x, r_)  return x:map(exp  , r_) end
function M.log   (x, r_)  return x:map(log  , r_) end
function M.sqrt  (x, r_)  return x:map(sqrt , r_) end
function M.proj  (x, r_)  return x:map(proj , r_) end
function M.sin   (x, r_)  return x:map(sin  , r_) end
function M.cos   (x, r_)  return x:map(cos  , r_) end
function M.tan   (x, r_)  return x:map(tan  , r_) end
function M.sinh  (x, r_)  return x:map(sinh , r_) end
function M.cosh  (x, r_)  return x:map(cosh , r_) end
function M.tanh  (x, r_)  return x:map(tanh , r_) end
function M.asin  (x, r_)  return x:map(asin , r_) end
function M.acos  (x, r_)  return x:map(acos , r_) end
function M.atan  (x, r_)  return x:map(atan , r_) end
function M.asinh (x, r_)  return x:map(asinh, r_) end
function M.acosh (x, r_)  return x:map(acosh, r_) end
function M.atanh (x, r_)  return x:map(atanh, r_) end
function M.unm   (x, r_)  return x:map(unm  , r_) end
function M.mod   (x, y, r_)  return x:map2(y, mod, r_) end -- TODO
function M.pow   (x, y, r_)  return x:map2(y, pow, r_) end -- TODO

function M.__eq (x, y)
  if iscalar(y) then
    for i=0,x:size()-1 do
      if x.data[i] ~= y then return false end
    end
    return true
  end

  if x:rows() ~= y:rows() or x:cols() ~= y:cols() then return false end
  for i=0,x:size()-1 do
    if x.data[i] ~= y.data[i] then return false end
  end
  return true
end

function M.add (x, y, r_)
  local r

  if isnum(x) then
    if ismat(y) then -- num + mat => vec + num
       r = r_ or matrix(y:sizes())
       assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
       clib.mad_vec_addn(y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num + cmat => cvec + num
       r = r_ or cmatrix(y:sizes())
       assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
       clib.mad_cvec_addn (y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    if isnum(y) then -- mat + num => vec + num
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_addn(x.data, y, r.data, r:size())
    elseif iscpx(y) then -- mat + cpx => vec + cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_addc(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then -- mat + mat => vec + vec
      r = r_ or matrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_add(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- mat + cmat => cvec + vec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_cvec_addv(y.data, x.data, r.data, r:size())
    else goto invalid end
    return r
  end

  if iscmat(x) then
    if isnum(y) then -- cmat + num => cvec + num
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_vec_addn(x.data, y, r.data, r:size())
    elseif iscpx(y) then -- cmat + cpx => cvec + cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_vec_addc(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then -- cmat + mat => cvec + vec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_vec_add(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- cmat + cmat => cvec + cvec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_addv(y.data, x.data, r.data, r:size())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (+) operands")
end

function M.sub (x, y, r_)
  local r

  if isnum(x) then
    if ismat(y) then -- num - mat => num - vec
      r = r_ or matrix(y:sizes())
      assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_subn (y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num - cmat => num - cvec
      r = r_ or cmatrix(y:sizes())
      assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_subn (y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    if isnum(y) then -- mat - num => vec + -num
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_addn(x.data, -y, r.data, r:size())
    elseif iscpx(y) then -- mat - cpx => vec + -cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_addc(x.data, -y.re, -y.im, r.data, r:size())
    elseif ismat(y) then -- mat - mat => vec - vec
      r = r_ or matrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_sub(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- mat - cmat => vec - cvec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_cvec_subv(x.data, y.data, r.data, r:size())
    else goto invalid end
    return r
  end

  if iscmat(x) then
    if isnum(y) then -- cmat - num => cvec + -num
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_addn(x.data, -y, r.data, r:size())
    elseif iscpx(y) then -- cmat - cpx => cvec + -cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_addc(x.data, -y.re, -y.im, r.data, r:size())
    elseif ismat(y) then -- cmat - mat => cvec - vec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_subv(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- cmat - cmat => cvec - cvec
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_sub(x.data, y.data, r.data, r:size())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (-) operands")
end

function M.mul (x, y, r_)
  local r

  if isnum(x) then
    if ismat(y) then -- num * mat => vec * num
      r = r_ or matrix(y:sizes())
      assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_muln (y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num * cmat => cvec * num
      r = r_ or cmatrix(y:sizes())
      assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_muln (y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    if isnum(y) then -- mat * num => vec * num
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_muln(x.data, y, r.data, r:size())
    elseif iscpx(y) then -- mat * cpx => vec * cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_mulc(x.data, y.re, y.im, r.data, r:size())
    elseif isvec(y) then -- mat * vec
      r = r_ or vector(x:rows())
      assert(x:cols() == y:size(),                   "incompatible matrix-vector sizes")
      assert(x:rows() == r:size(),                   "incompatible matrix-vector sizes")
      clib.mad_mat_muln(x.data, y.data, r.data, x:rows(), x:cols())
    elseif iscvec(y) then -- mat * cvec
      r = r_ or cvector(x:rows())
      assert(x:cols() == y:size(),                   "incompatible matrix-cvector sizes")
      assert(x:rows() == r:size(),                   "incompatible matrix-cvector sizes")
      clib.mad_mat_mulc(x.data, y.data, r.data, x:rows(), x:cols())
    elseif ismat(y) then -- mat * mat
      r = r_ or matrix(x:rows(), y:cols())
      assert(x:cols() == y:rows(),                          "incompatible matrix sizes")
      assert(x:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_mul(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    elseif iscmat(y) then -- mat * cmat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows(),                          "incompatible matrix sizes")
      assert(x:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_mulm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    else goto invalid end
   return r
  end

  if iscmat(x) then
    if isnum(y) then -- cmat * num => cvec * num
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_muln(x.data, y, r.data, r:size())
    elseif iscpx(y) then -- cmat * cpx => cvec * cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_mulc(x.data, y.re, y.im, r.data, r:size())
    elseif isvec(y) then -- cmat * vec
      r = r_ or cvector(x:rows())
      assert(x:cols() == y:size(),                  "incompatible cmatrix-vector sizes")
      assert(x:rows() == r:size(),                  "incompatible cmatrix-vector sizes")
      clib.mad_cmat_muln(x.data, y.data, r.data, x:rows(), x:cols())
    elseif iscvec(y) then -- cmat * cvec
      r = r_ or cvector(x:rows())
      assert(x:cols() == y:size(),                  "incompatible cmatrix-cvector sizes")
      assert(x:rows() == r:size(),                  "incompatible cmatrix-cvector sizes")
      clib.mad_cmat_mulc(x.data, y.data, r.data, x:rows(), x:cols())
    elseif ismat(y) then -- cmat * mat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows(),                          "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_mulm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    elseif iscmat(y) then -- cmat * cmat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows(),                          "incompatible cmatrix sizes")
      assert(x:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_mul(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (*) operands")
end

local lapack = require 'lapack'
local T = ffi.new('char[1]', 'T')
local Ms = ffi.new('int32_t[1]')
local Ns = ffi.new('int32_t[1]')
local Ks = ffi.new('int32_t[1]')
local ZERO = ffi.new('double[1]')
local ONE  = ffi.new('double[1]', 1.0)

function M.div (x, y, r_)
  local r

  if isnum(x) then
    if ismat(y) then -- num / mat => num * inv(mat)
      error("num/mat: NYI matrix inverse")
    elseif ismat(y) then -- num / cmat => num * inv(cmat)
      error("num/cmat: NYI cmatrix inverse")
    else goto invalid end
  end

  if ismat(x) then
    if isnum(y) then -- mat / num => vec * (1/num)
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_muln(x.data, 1/y, r.data, r:size())
    elseif iscpx(y) then -- mat / cpx => vec * (1/cpx)
      r, y = r_ or cmatrix(x:sizes()), 1/y
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_mulc(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then -- mat / mat => mat * inv(mat)
--      Ms[0], Ns[0], Ks[0] = r:rows(), r:cols(), x:cols()
--      lapack.dgemm_(T, T, Ms, Ns, Ks, ONE, x.data, Ns, y.data, Ns, ZERO, r.data, Ns)
      error("mat/mat: NYI matrix inverse")
    elseif iscmat(y) then -- mat / cmat => mat * inv(cmat)
      error("mat/cmat: NYI matrix inverse")
    else goto invalid end
    return r
  end

  if iscmat(x) then
    if isnum(y) then -- cmat / num => cvec * (1/num)
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_muln(x.data, 1/y, r.data, r:size())
    elseif iscpx(y) then -- cmat / cpx => cvec * (1/cpx)
      r, y = r_ or cmatrix(x:sizes()), 1/y
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_mulc(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then -- cmat / mat => cmat * inv(mat)
      error("cmat/mat: NYI matrix inverse")
    elseif iscmat(y) then -- cmat / cmat => cmat * inv(cmat)
      error("cmat/cmat: NYI cmatrix inverse")
    else goto invalid end
    return r
  end

::invalid:: error("incompatible matrix (/) operands")
end

function M.schur_prod (x, y, r_)
  local nr, nc, r = x:sizes()

  if ismat(x) then
    if ismat(y) then -- mat .* mat => vec * vec
      r = r_ or matrix(nr, nc)
      assert(nr == y:rows() and nc == y:cols(), "incompatible matrix sizes")
      assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_mul(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- mat .* cmat => cvec * vec
      r = r_ or cmatrix(nr, nc)
      assert(nr == y:rows() and nc == y:cols(), "incompatible matrix sizes")
      assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
      clib.mad_cvec_mulv(y.data, x.data, r.data, r:size())
    else goto invalid end
    return r
  end

  if iscmat(x) then
    if ismat(y) then -- mat .* cmat => cvec * vec
      r = r_ or cmatrix(nr, nc)
      assert(nr == y:rows() and nc == y:cols(), "incompatible matrix sizes")
      assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
      clib.mad_cvec_mulv(y.data, x.data, r.data, nr*nc)
    elseif iscmat(y) then -- cmat .* cmat => cvec * cvec
      r = r_ or cmatrix(nr, nc)
      assert(nr == y:rows() and nc == y:cols(), "incompatible matrix sizes")
      assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
      clib.mad_cvec_mul(x.data, y.data, r.data, nr*nc)
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (.*) operands")
end

function M.tostring (x, sep_, lsep_)
  local nr, nc = x:sizes()
  local r, c = {}, {}
  for i=1,nr do
    for j=1,nc do c[j] = tostring(x:get(i,j)) end
    r[i] = table.concat(c, sep_ or ' ')
  end
  return table.concat(r, lsep_ or '\n')
end

local tbl_new = require 'table.new'

function M.totable(x, r_)
  local nr, nc = x:sizes()
  local r = r_ or tbl_new(nr,0)
  assert(type(r) == 'table', "invalid argument, table expected")
  for i=1,nr do
    local c = r[i] or tbl_new(nc,0)
    assert(type(c) == 'table', "invalid argument, table of tables expected")
    for j=1,nc do c[j] = x:get(i,j) end
    r[i] = c
  end
  return r
end

M.__unm      = M.unm
M.__add      = M.add
M.__sub      = M.sub
M.__mul      = M.mul
M.__div      = M.div
M.__mod      = M.mod
M.__pow      = M.pow
M.__tostring = M.tostring
M.__index    = M

ffi.metatype( 'matrix_t', M)
ffi.metatype('cmatrix_t', M)

-- END -------------------------------------------------------------------------
return matrix
