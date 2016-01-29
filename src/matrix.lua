--[=[
 o----------------------------------------------------------------------------o
 |
 | Matrix module (real)
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
  
  Purpose:
  - provides full set of functions and operations on real and complex matrices

  Information:
  - real and complex matrices are implemented by this module

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  matrix

SYNOPSIS
  local matrix = require 'matrix'
  local m1 = matrix(3)                      -- column matrix = matrix(3,1)
  local m2 = matrix(2,3)
  local m3 = matrix {{1,2},{3,4},{5,6}}
  local m4 = matrix {1,2,3,4,5,6}           -- column matrix = {{1},{2},...}
  local m5 = matrix {{1,2,3,4,5,6}}         -- row matrix
  local m6 = m1:transpose()                 -- row matrix
  local I6 = matrix(6):ones()               -- 6x6 identity

DESCRIPTION
  The module matrix implements the operators and math functions on matrices:
    (minus) -, +, -, *, /, %, ^, ==, #, [], ..,
    unm, add, sub, mul, div, mod, pow, emul, ediv,
    rows, cols, size, sizes, get, set, get0, set0,
    zeros, ones, eye, unit, fill, copy,
    get_row, get_col, get_diag, get_sub,
    set_row, set_col, set_diag, set_sub,
    transpose, t, trans, ctrans,
    real, imag, conj, norm, angle, trace, tr,
    dot, inner, cross, mixed, outer,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, foreach, map, map2, maps,
    solve, svd, eigen,
    concat, reshape, tostring, totable, fromtable,
    check_bounds.

REMARK:
  By default, check_bounds is true.

RELATIONS (3D GEOMETRY)
  inner prod:   u'.v = |u|.|v| cos(u^v)
  cross prod:   uxv = |u|.|v| sin(u^v) \vec{n}
  mixed prod:   (uxv)'.w = u'.(vxw) = det(u,v,w)
  outer prod:   u.v' = matrix
  dble xprod:   ux(vxw) = (u.w) \vec{v} - (u.v) \vec{w}
                (uxv)xw = (u.w) \vec{v} - (v.w) \vec{u}
  norm      :   |u| = sqrt(u'.u)
  angle     :   u^v = acos(u'.v / |u|.|v|)  in [0,pi] (or [-pi,pi] if n)
  unit      :   u / |u|
  projection:   u'.v
  projector :   I -   u.u' / u'.u
  reflector :   I - 2 u.u' / u'.u
  area      :   |uxv|
  volume    :   |(uxv)'.w|
  unitary   :   |u| = 1
  orthogonal:   u'.v = 0
  collinear :   |uxv| = 0
  coplanar  :   (uxv)'.w = 0

RETURN VALUES
  The constructor of matrices

SEE ALSO
  gmath, complex, matrix, cmatrix
]]
 
-- modules -------------------------------------------------------------------o

local ffi     = require 'ffi'
local clib    = require 'cmad'
local gmath   = require 'gmath'
local xmatrix = require 'xmatrix'

local tbl_new = require 'table.new'

-- locals --------------------------------------------------------------------o

local isnum, iscpx, iscal,
      ismat, iscmat, isamat,
      real, imag, conj, ident, min,
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, mod, pow, tostring = 
      gmath.is_number, gmath.is_complex, gmath.is_scalar,
      gmath.is_matrix, gmath.is_cmatrix, gmath.isa_matrix,
      gmath.real, gmath.imag, gmath.conj, gmath.ident, gmath.min,
      gmath.abs, gmath.arg, gmath.exp, gmath.log, gmath.sqrt, gmath.proj,
      gmath.sin, gmath.cos, gmath.tan, gmath.sinh, gmath.cosh, gmath.tanh,
      gmath.asin, gmath.acos, gmath.atan, gmath.asinh, gmath.acosh, gmath.atanh,
      gmath.unm, gmath.mod, gmath.pow, gmath.tostring

local istype, cast, sizeof, fill = ffi.istype, ffi.cast, ffi.sizeof, ffi.fill

local cres = ffi.new 'complex[1]'

-- FFI type constructors
local matrix  = xmatrix.matrix
local cmatrix = xmatrix.cmatrix

-- implementation ------------------------------------------------------------o

local function unsafe_get0 (x, i, j) return x.data[ i    * x.nc +  j   ]     end
local function unsafe_set0 (x, i, j, e)     x.data[ i    * x.nc +  j   ] = e end
local function unsafe_get  (x, i, j) return x.data[(i-1) * x.nc + (j-1)]     end
local function unsafe_set  (x, i, j, e)     x.data[(i-1) * x.nc + (j-1)] = e end

local function safe_get0 (x, i, j)
  assert(0 <= i and i < x.nr and 0 <= j and j < x.nc, "0-index out of bounds")
  return unsafe_get0(x, i, j)
end

local function safe_set0 (x, i, j, e)
  assert(0 <= i and i < x.nr and 0 <= j and j < x.nc, "0-index out of bounds")
  unsafe_set0(x, i, j, e)
end

local function safe_get (x, i, j)   
  assert(1 <= i and i <= x.nr and 1 <= j and j <= x.nc, "1-index out of bounds")
  return unsafe_get(x, i, j)
end

local function safe_set (x, i, j, e)
  assert(1 <= i and i <= x.nr and 1 <= j and j <= x.nc, "1-index out of bounds")
  unsafe_set(x, i, j, e)
end

function M.check_bounds (self, flag)
  local old = M.get0 == safe_get0
  if flag == true then
    M.get0, M.get = safe_get0, safe_get
    M.set0, M.set = safe_set0, safe_set
  else
    M.get0, M.get = unsafe_get0, unsafe_get
    M.set0, M.set = unsafe_set0, unsafe_set
  end
  return old;
end

M.check_bounds(nil, true) -- default

function M.rows   (x) return x.nr        end
function M.cols   (x) return x.nc        end
function M.size   (x) return x.nr * x.nc end
function M.sizes  (x) return x.nr , x.nc end
function M.tsizes (x) return x.nc , x.nr end

function M.reshape (x, nr, nc)
  assert(nr*nc <= x:size(), "incompatible matrix sizes")
  x.nr, x.nc = nr, nc
  return x
end

function M.zeros (x)
  fill(x.data, sizeof(ismat(x) and 'double' or 'complex', x:size()))
  return x
end

function M.ones (x, e_)
  local n, e = x:size(), e_ or 1
  for i=0,n-1 do x.data[i] = e end
  return x
end

function M.eye (x, e_) -- zeros + diagonal
  x:zeros()
  local n, e = min(x:sizes()), e_ or 1
  for i=0,n-1 do x:set0(i,i, e) end
  return x
end

function M.unit(x)
  local n = x:norm()
  assert(n ~= 0, "null matrix")
  if n ~= 1 then x:div(n, x) end
  return x
end

function M.get_col (x, jc, r_)
  local j, nr, nc = jc-1, x:sizes()
  assert(0 <= j and j < nc, "column index out of bounds")
  local r = r_ or ismat(x) and matrix(nr,1) or cmatrix(nr,1)
  assert(not isamat(r) or nr == r:size(), "incompatible matrix sizes")
  for i=0,nr-1 do r[i+1] = x:get0(i,j) end
  return r
end

function M.get_row (x, ir, r_)
  local i, nr, nc = ir-1, x:sizes()
  assert(0 <= i and i < nr, "row index out of bounds")
  local r = r_ or ismat(x) and matrix(1,nc) or cmatrix(1,nc)
  assert(not isamat(r) or nc == r:size(), "incompatible matrix sizes")
  for j=0,nc-1 do r[j+1] = x:get0(i,j) end
  return r
end

function M.get_diag (x, r_)
  local n = min(x:sizes())
  local r = r_ or ismat(x) and matrix(n,1) or cmatrix(n,1)
  assert(not isamat(r) or n == r:size(), "incompatible matrix sizes")
  for i=0,n-1 do r[i+1] = x:get0(i,i) end
  return r
end

function M.get_sub (x, rg, r_)
  local i, j, nr, nc = rg[1]-1, rg[2]-1, rg[3], rg[4]
  assert(0 <= i and i+nr <= x:rows(), "row index out of bounds")
  assert(0 <= j and j+nc <= x:cols(), "column index out of bounds")
  local r = r_ or ismat(x) and matrix(nr,nc) or cmatrix(nr,nc)
  assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
  for ii=0,nr-1 do
  for jj=0,nc-1 do r:set0(ii,jj, x:get0(i+ii,j+jj)) end end
  return r
end

function M.set_col (x, jc, y)
  local j, nr, nc = jc-1, x:sizes() 
  assert(0 <= j and j < nc, "column index out of bounds")
  local n = min(nr, #y)
  for i=0,n-1 do x:set0(i,j, y[i+1]) end
  return x
end

function M.set_row (x, ir, y)
  local i, nr, nc = ir-1, x:sizes()
  assert(0 <= i and i < nr, "row index out of bounds")
  local n = min(nc, #y)
  for j=0,n-1 do x:set0(i,j, y[j+1]) end
  return x
end

function M.set_diag (x, y)
  local n = min(#y, x:sizes())
  for i=0,n-1 do x:set0(i,i, y[i+1]) end
  return x
end

function M.set_sub (x, rg, y)
  local i, j, nr, nc = rg[1]-1, rg[2]-1, y:sizes()
  assert(0 <= i and i+nr <= x:rows(), "row index out of bounds")
  assert(0 <= j and j+nc <= x:cols(), "column index out of bounds")
  for ii=0,nr-1 do
  for jj=0,nc-1 do x:set0(i+ii,j+jj, y:get0(ii,jj)) end end
  return r
end

function M.foldl (x, r, f)
  local nr, nc = x:sizes()
  if nr == r:rows() then
    for i=0,nr-1 do
      local ri = r:get0(i, 0)
      for j=0,nc-1 do ri = f(ri, x:get0(i,j)) end
      r:set0(i, 0, ri)
    end
  elseif nc == r:cols() then
    for j=0,nc-1 do
      local rj = r:get0(0, j)
      for i=0,nr-1 do rj = f(rj, x:get0(i,j)) end
      r:set0(0, j, rj)
    end
  else
    error("incompatible matrix sizes")
  end
  return r
end

function M.foldr (x, r, f)
  local nr, nc = x:sizes()
  if nr == r:rows() then
    for i=0,nr-1 do
      local ri = r:get0(i, 0)
      for j=0,nc-1 do ri = f(x:get0(i,j), ri) end
      r:set0(i, 0, ri)
    end
  elseif nc == r:cols() then
    for j=0,nc-1 do
      local rj = r:get0(0, j)
      for i=0,nr-1 do rj = f(x:get0(i,j), rj) end
      r:set0(0, j, rj)
    end
  else
    error("incompatible matrix sizes")
  end
  return r
end

function M.foreach (x, f)
  for i=0,x:size()-1 do f(x.data[i]) end
  return x
end

function M.map (x, f, r_)
  local r0 = f(x.data[0])
  local r = r_ or isnum(r0) and matrix(x:sizes()) or cmatrix(x:sizes())
  assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
  r.data[0] = r0
  for i=1,x:size()-1 do r.data[i] = f(x.data[i]) end
  return r
end

function M.map2 (x, y, f, r_)
  assert(isamat(y), "matrix expected for 2nd argument")
  local r0 = f(x.data[0], y.data[0])
  local r = r_ or isnum(r0) and matrix(x:sizes()) or cmatrix(x:sizes())
  assert(x:rows() == y:rows() and x:cols() == y:cols() and
         x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
  r.data[0] = r0
  for i=1,x:size()-1 do r.data[i] = f(x.data[i], y.data[i]) end
  return r
end

function M.maps (x, y, f, r_)
  if ismat(y) then return x:map2(y,f,r_) end
  assert(iscal(y), "scalar expected for 2nd argument")
  local r0 = f(x.data[0], y)
  local r = r_ or isnum(r0) and matrix(x:sizes()) or cmatrix(x:sizes())
  assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
  r.data[0] = r0
  for i=1,x:size()-1 do r.data[i] = f(x.data[i], y) end
  return r
end

function M.conjugate (x, r_)
  local nr, nc = x:sizes()
  local r = r_ or ismat(x) and matrix(nc, nr) or cmatrix(nc, nr)

  assert(nr == r:cols() and nc == r:rows(), "incompatible matrix sizes")
  if ismat(x) and x.data ~= r.data then
    clib.mad_vec_copy(x.data, r.data, nr*nc)  -- copy
  elseif iscmat(x) then
    clib.mad_cvec_conj(x.data, r.data, nr*nc) -- conjugate
  end
  return r
end

M.conj = M.conjugate -- shortcut

function M.transpose (x, r_, c_)
  local nr, nc = x:sizes()
  local r = r_ or ismat(x) and matrix(nc, nr) or cmatrix(nc, nr)
  local c = c_ or true
  assert(nr == r:cols() and nc == r:rows(), "incompatible matrix sizes")
  if ismat(x) then
    clib.mad_mat_trans(x.data, r.data, nr, nc)   -- transpose
  elseif c == true then
    clib.mad_cmat_ctrans(x.data, r.data, nr, nc) -- conjugate transpose
  else
    clib.mad_cmat_trans(x.data, r.data, nr, nc)  -- transpose (no conjugate)
  end
  return r
end

M.t = M.transpose -- shortcut
M.trans  = function (x, r_) return x:t(r_, false) end -- never  conjugate
M.ctrans = function (x, r_) return x:t(r_, true ) end -- always conjugate

function M.trace (x)
  local n, r = min(x:sizes()), 0
  for i=0,n-1 do r = r + x:get0(i,i) end
  return r
end

M.tr = M.trace -- shortcut

function M.inner (x, y)
  -- (x:t() * y):tr() without temporaries
  assert(isamat(y), "matrix expected for 2nd argument")
  assert(x:rows() == y:rows(), "incompatible matrix sizes")

  if ismat(x) then
    if ismat(y) then
      return clib.mad_mat_dot(x.data, y.data, x:cols(), y:cols(), x:rows())
    else
      clib.mad_mat_dotm_r(x.data, y.data, cres, x:cols(), y:cols(), x:rows())
      return cres[0]
    end
  end

  if iscmat(x) then
    if ismat(y) then
      clib.mad_cmat_dotm_r(x.data, y.data, cres, x:cols(), y:cols(), x:rows())
    else
      clib.mad_cmat_dot_r(x.data, y.data, cres, x:cols(), y:cols(), x:rows())
    end
    return cres[0]
  end

  error("invalid matrix dot operands")
end

M.dot = M.inner -- shortcut

function M.cross (x, y, r_)
  assert(isamat(y), "matrix expected for 2nd argument")
  local nr, nc = x:sizes()
  assert(nr == 3, "invalid matrix sizes")
  local r = r_ or ismat(x) and ismat(y) and matrix(3,nc) or cmatrix(3,nc)
  assert(nr == y:rows() and nc == y:cols() and
         nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
  local i0, i1, i2 = 0, nc, 2*nc
  for i=0,nc-1 do
    r.data[i0] = x.data[i1] * y.data[i2] - x.data[i2] * y.data[i1]
    r.data[i1] = x.data[i2] * y.data[i0] - x.data[i0] * y.data[i2]
    r.data[i2] = x.data[i0] * y.data[i1] - x.data[i1] * y.data[i0]
    i0, i1, i2 = i0+1, i1+1, i2+1
  end
  return r
end

function M.mixed (x, y, z, r_)
  -- x:cross(y):dot(z) without temporary
  assert(isamat(y) and isamat(z), "matrices expected for 2nd and 3rd arguments")
  local nr, nc = x:sizes()
  assert(nr == 3, "invalid matrix sizes")
  local r
  if nc == 1 then -- single mixed product (r_ ignored)
    assert(nr == y:rows() and 1 == y:cols() and
           nr == z:rows() and 1 == z:cols(), "incompatible matrix sizes")
    r = conj(x.data[1] * y.data[2] - x.data[2] * y.data[1]) * z.data[0] +
        conj(x.data[2] * y.data[0] - x.data[0] * y.data[2]) * z.data[1] +
        conj(x.data[0] * y.data[1] - x.data[1] * y.data[0]) * z.data[2]
  else -- multiple mixed products
    r = r_ or ismat(x) and ismat(y) and matrix(nc,1) or cmatrix(nc,1)
    assert(nr == y:rows() and nc == y:cols() and
           nr == z:rows() and nc == z:cols() and
           nc == r:rows() and 1  == r:cols(), "incompatible matrix sizes")
    local i0, i1, i2 = 0, nc, 2*nc
    for i=0,nc-1 do
      r.data[i] = conj(x.data[i1] * y.data[i2] - x.data[i2] * y.data[i1]) * z.data[i0] +
                  conj(x.data[i2] * y.data[i0] - x.data[i0] * y.data[i2]) * z.data[i1] +
                  conj(x.data[i0] * y.data[i1] - x.data[i1] * y.data[i0]) * z.data[i2]
      i0, i1, i2 = i0+1, i1+1, i2+1
    end
  end
  return r
end

function M.outer (x, y, r_)
  -- x * y:t() without temporary
  assert(isamat(y), "matrix expected for 2nd argument")
  local nr, nc = x:rows(), y:rows()
  assert(x:cols() == 1 and y:cols() == 1, "invalid matrix sizes")
  local r = r_ or ismat(x) and ismat(y) and matrix(nr,nc) or cmatrix(nr,nc)
  assert(nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
  for i=0,nr-1 do
    for j=0,nc-1 do r:set0(i,j, x.data[i] * conj(y.data[j])) end
  end
  return r
end

function M.norm (x)
  -- Frobenius norm (consistent with inner product)
  if ismat(x) then
    return sqrt(clib.mad_vec_dot(x.data, x.data, x:size()))
  else
    clib.mad_cvec_dot_r(x.data, x.data, cres, x:size())
    return sqrt(cres[0])
  end
end

function M.angle (x, y, n_)
  local w = x:inner(y)
  local v = x:norm() * y:norm()
  assert(v ~= 0, "null vector") -- convention: return pi/2 ?
  local a = acos(w / v) -- [0, pi]
  if n_ and x:mixed(y, n_) < 0 then a = -a end -- [-pi, pi]
  return a
end

function M.fill  (x, e )    return x:map (function() return e end, x) end
function M.copy  (x, r_)    return x:map (ident , r_) end
function M.real  (x, r_)    return x:map (real  , r_) end
function M.imag  (x, r_)    return x:map (imag  , r_) end
function M.unm   (x, r_)    return x:map (unm   , r_) end
function M.abs   (x, r_)    return x:map (abs   , r_) end
function M.arg   (x, r_)    return x:map (arg   , r_) end
function M.exp   (x, r_)    return x:map (exp   , r_) end
function M.log   (x, r_)    return x:map (log   , r_) end
function M.sqrt  (x, r_)    return x:map (sqrt  , r_) end
function M.proj  (x, r_)    return x:map (proj  , r_) end
function M.sin   (x, r_)    return x:map (sin   , r_) end
function M.cos   (x, r_)    return x:map (cos   , r_) end
function M.tan   (x, r_)    return x:map (tan   , r_) end
function M.sinh  (x, r_)    return x:map (sinh  , r_) end
function M.cosh  (x, r_)    return x:map (cosh  , r_) end
function M.tanh  (x, r_)    return x:map (tanh  , r_) end
function M.asin  (x, r_)    return x:map (asin  , r_) end
function M.acos  (x, r_)    return x:map (acos  , r_) end
function M.atan  (x, r_)    return x:map (atan  , r_) end
function M.asinh (x, r_)    return x:map (asinh , r_) end
function M.acosh (x, r_)    return x:map (acosh , r_) end
function M.atanh (x, r_)    return x:map (atanh , r_) end
function M.mod   (x, y, r_) return x:maps(y, mod, r_) end
function M.pow   (x, y, r_) return x:maps(y, pow, r_) end

--[[ TODO
  det
  mexp, mlog, msqrt, mpow
  __eq with tol
]]

function M.__eq (x, y)
  if iscal(y) then
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
    r = r_ or iscmat(y) and cmatrix(y:sizes()) or matrix(y:sizes())
    assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
    if ismat(y) then      -- num + mat => vec + num
      clib.mad_vec_addn(y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num + cmat => cvec + num
      clib.mad_cvec_addn(y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    r = r_ or (iscpx(y) or iscmat(y)) and cmatrix(x:sizes()) or matrix(x:sizes())
    assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
    if isnum(y) then      -- mat + num => vec + num
      clib.mad_vec_addn(x.data, y, r.data, r:size())
    elseif iscpx(y) then  -- mat + cpx => vec + cpx
      clib.mad_vec_addc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- mat + mat => vec + vec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      clib.mad_vec_add(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- mat + cmat => cvec + vec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      clib.mad_cvec_addv(y.data, x.data, r.data, r:size())
    else goto invalid end
    return r
  end

  if iscmat(x) then
    r = r_ or cmatrix(x:sizes())
    assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
    if isnum(y) then      -- cmat + num => cvec + num
      clib.mad_vec_addn(x.data, y, r.data, r:size())
    elseif iscpx(y) then  -- cmat + cpx => cvec + cpx
      clib.mad_vec_addc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- cmat + mat => cvec + vec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      clib.mad_vec_addv(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- cmat + cmat => cvec + cvec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_add(y.data, x.data, r.data, r:size())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (+) operands")
end

function M.sub (x, y, r_)
  local r

  if isnum(x) then
    r = r_ or (iscpx(y) or iscmat(y)) and cmatrix(y:sizes()) or matrix(y:sizes())
    assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
    if ismat(y) then      -- num - mat => num - vec
      clib.mad_vec_subn(y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num - cmat => num - cvec
      clib.mad_cvec_subn(y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    r = r_ or (iscpx(y) or iscmat(y)) and cmatrix(x:sizes()) or matrix(x:sizes())
    assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
    if isnum(y) then      -- mat - num => vec + -num
      clib.mad_vec_addn(x.data, -y, r.data, r:size())
    elseif iscpx(y) then  -- mat - cpx => vec + -cpx
      clib.mad_vec_addc_r(x.data, -y.re, -y.im, r.data, r:size())
    elseif ismat(y) then  -- mat - mat => vec - vec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      clib.mad_vec_sub(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- mat - cmat => vec - cvec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible matrix sizes")
      clib.mad_cvec_subv(x.data, y.data, r.data, r:size())
    else goto invalid end
    return r
  end

  if iscmat(x) then
    r = r_ or cmatrix(x:sizes())
    assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
    if isnum(y) then      -- cmat - num => cvec + -num
      clib.mad_cvec_addn(x.data, -y, r.data, r:size())
    elseif iscpx(y) then  -- cmat - cpx => cvec + -cpx
      clib.mad_cvec_addc_r(x.data, -y.re, -y.im, r.data, r:size())
    elseif ismat(y) then  -- cmat - mat => cvec - vec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_subv(x.data, y.data, r.data, r:size())
    elseif iscmat(y) then -- cmat - cmat => cvec - cvec
      assert(x:rows() == y:rows() and x:cols() == y:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_sub(x.data, y.data, r.data, r:size())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (-) operands")
end

function M.mul (x, y, r_)
  local r

  if isnum(x) then
    r = r_ or (iscpx(y) or iscmat(y)) and cmatrix(y:sizes()) or matrix(y:sizes())
    assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
    if ismat(y) then      -- num * mat => vec * num
      clib.mad_vec_muln(y.data, x, r.data, r:size())
    elseif iscmat(y) then -- num * cmat => cvec * num
      clib.mad_cvec_muln(y.data, x, r.data, r:size())
    else goto invalid end
    return r
  end

  if ismat(x) then
    if isnum(y) then      -- mat * num => vec * num
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_muln(x.data, y, r.data, r:size())
    elseif iscpx(y) then  -- mat * cpx => vec * cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_mulc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- mat * mat
      r = r_ or matrix(x:rows(), y:cols())
      assert(x:cols() == y:rows() and x:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_mul(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    elseif iscmat(y) then -- mat * cmat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows() and x:rows() == r:rows() and y:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_mulm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    else goto invalid end
   return r
  end

  if iscmat(x) then
    if isnum(y) then      -- cmat * num => cvec * num
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_muln(x.data, y, r.data, r:size())
    elseif iscpx(y) then  -- cmat * cpx => cvec * cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_mulc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- cmat * mat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows() and x:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_mulm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    elseif iscmat(y) then -- cmat * cmat
      r = r_ or cmatrix(x:rows(), y:cols())
      assert(x:cols() == y:rows() and x:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_mul(x.data, y.data, r.data, r:rows(), r:cols(), x:cols())
    else goto invalid end
    return r
  end

::invalid:: error("invalid matrix (*) operands")
end

function M.div (x, y, r_, rcond_)
  local r

  if isnum(x) then
    r = r_ or iscmat(y) and cmatrix(y:tsizes()) or matrix(y:tsizes())
    assert(y:rows() == r:cols() and y:cols() == r:rows(), "incompatible matrix sizes")
    if ismat(y) then      -- num / mat
      clib.mad_mat_invn (y.data, x, r.data, y:rows(), y:cols(), rcond_ or -1)
    elseif iscmat(y) then -- num / cmat
      clib.mad_cmat_invn(y.data, x, r.data, y:rows(), y:cols(), rcond_ or -1)
    else goto invalid end
    return r
  end

  if ismat(x) then
    if isnum(y) then       -- mat / num => vec / num
      r = r_ or matrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_divn(x.data, y, r.data, r:size())
    elseif iscpx(y) then   -- mat / cpx => vec / cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")
      clib.mad_vec_divc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- mat / mat => X * Y^-1 => [m x n] = [m x p] * [p x n] => X:[m x p], Y:[n x p]
      r = r_ or matrix(x:rows(), y:rows())
      assert(x:cols() == y:cols() and x:rows() == r:rows() and y:rows() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_div(x.data, y.data, r.data, r:rows(), r:cols(), x:cols(), rcond_ or -1)
    elseif iscmat(y) then  -- mat / cmat
      r = r_ or cmatrix(x:rows(), y:rows())
      assert(x:cols() == y:cols() and x:rows() == r:rows() and y:rows() == r:cols(), "incompatible matrix sizes")
      clib.mad_mat_divm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols(), rcond_ or -1)
    else goto invalid end
    return r
  end

  if iscmat(x) then
    if isnum(y) then       -- cmat / num => cvec / num
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_muln(x.data, 1/y, r.data, r:size())
    elseif iscpx(y) then   -- cmat / cpx => cvec / cpx
      r = r_ or cmatrix(x:sizes())
      assert(x:rows() == r:rows() and x:cols() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cvec_divc_r(x.data, y.re, y.im, r.data, r:size())
    elseif ismat(y) then  -- cmat / mat
      r = r_ or cmatrix(x:rows(), y:rows())
      assert(x:cols() == y:cols() and x:rows() == r:rows() and y:rows() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_divm(x.data, y.data, r.data, r:rows(), r:cols(), x:cols(), rcond_ or -1)
    elseif iscmat(y) then  -- cmat / cmat
      r = r_ or cmatrix(x:rows(), y:rows())
      assert(x:cols() == y:cols() and x:rows() == r:rows() and y:rows() == r:cols(), "incompatible cmatrix sizes")
      clib.mad_cmat_div(x.data, y.data, r.data, r:rows(), r:cols(), x:cols(), rcond_ or -1)
    else goto invalid end
    return r
  end

::invalid:: error("incompatible matrix (/) operands")
end

function M.emul (x, y, r_)
  assert(isamat(y), "matrix expected for 2nd argument")
  local r = r_ or ismat(x) and ismat(y) and matrix(x:sizes()) or cmatrix(x:sizes())
  assert(x:rows() == y:rows() and x:cols() == y:cols() and
         x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")

  if ismat(x) then
    if ismat(y) then    -- mat .* mat => vec * vec
      clib.mad_vec_mul(x.data, y.data, r.data, r:size())
    else                -- mat .* cmat => cvec * vec
      clib.mad_cvec_mulv(y.data, x.data, r.data, r:size())
    end
    return r
  end

  if iscmat(x) then
    if ismat(y) then    -- mat .* cmat => cvec * vec
      clib.mad_cvec_mulv(y.data, x.data, r.data, r:size())
    else                -- cmat .* cmat => cvec * cvec
      clib.mad_cvec_mul(x.data, y.data, r.data, r:size())
    end
    return r
  end

  error("invalid matrix (.*) operands")
end

function M.ediv (x, y, r_)
  assert(isamat(y), "matrix expected for 2nd argument")
  local r = r_ or ismat(x) and ismat(y) and matrix(x:sizes()) or cmatrix(x:sizes())
  assert(x:rows() == y:rows() and x:cols() == y:cols() and
         x:rows() == r:rows() and x:cols() == r:cols(), "incompatible matrix sizes")

  if ismat(x) then
    if ismat(y) then    -- mat ./ mat => vec / vec
      clib.mad_vec_div(x.data, y.data, r.data, r:size())
    else                -- mat ./ cmat => vec / cvec
      clib.mad_vec_divv(x.data, y.data, r.data, r:size())
    end
    return r
  end

  if iscmat(x) then
    if ismat(y) then    -- mat ./ cmat => cvec / vec
      clib.mad_cvec_divv(x.data, y.data, r.data, r:size())
    else                -- cmat ./ cmat => cvec / cvec
      clib.mad_cvec_div(x.data, y.data, r.data, r:size())
    end
    return r
  end

  error("invalid matrix (./) operands")
end

function M.svd (x)
  local nr, nc = x:sizes()

  if ismat(x) then
    local ru, rv, rs = matrix(nr,nr), matrix(nc,nc), matrix(min(nr,nc),1)
    local info = clib.mad_mat_svd(x.data, ru.data, rs.data, rv.data, x:rows(), x:cols())
    return ru, rs, rv, info
  else
    local ru, rv, rs = cmatrix(nr,nr), cmatrix(nc,nc), matrix(min(nr,nc),1)
    local info = clib.mad_cmat_svd(x.data, ru.data, rs.data, rv.data, x:rows(), x:cols())
    return ru, rs, rv, info
  end
end

function M.eigen (x)
  local nr, nc = x:sizes()
  assert(nr == nc, "matrix must be square")

  if ismat(x) then
    local wr, wi = matrix(nr, 1), matrix(nr, 1)
    local vl, vr = matrix(nr,nr), matrix(nr,nr)
    local info = clib.mad_mat_eigen(x.data, wr.data, wi.data, vl.data, vr.data, x:rows())
    return wr, wi, vl, vr, info
  else
    local w      = cmatrix(nr, 1)
    local vl, vr = cmatrix(nr,nr), cmatrix(nr,nr)
    local info = clib.mad_cmat_eigen(x.data, w.data, vl.data, vr.data, x:rows())
    return w, vl, vr, info
  end
end

function M.concat (x, y, v_, r_)
  assert(isamat(y), "matrix expected for 2nd argument")
  local nrx, ncx = x:sizes()
  local nry, ncy = y:sizes()
  if v_ then -- concat columns (vectical)
    local nr, nc = nrx + nry, ncx
    local r = r_ or ismat(x) and ismat(y) and matrix(nr,nc) or cmatrix(nr,nc)
    assert(ncx == ncy and nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
    local nx, ny = nrx * nc, nry * nc
    for i=0,nx-1 do r.data[i   ] = x.data[i] end
    for i=0,ny-1 do r.data[i+nx] = y.data[i] end
  else -- concat rows (horizontal)
    local nr, nc = nrx, ncx + ncy
    local r = r_ or ismat(x) and ismat(y) and matrix(nr,nc) or cmatrix(nr,nc)
    assert(nrx == nry and nr == r:rows() and nc == r:cols(), "incompatible matrix sizes")
    for i=0,nr do
      for j=0,ncx-1 do r:set0(i,j    , x:get0(i,j)) end
      for j=0,ncy-1 do r:set0(i,j+ncx, y:get0(i,j)) end
    end
  end
  return r
end

function M.tostring (x, sep_, lsep_)
  local nr, nc = x:sizes()
  local r, c = tbl_new(nr,0), tbl_new(nc,0)
  for i=0,nr-1 do
    for j=0,nc-1 do c[j+1] = tostring(x:get0(i,j)) end
    r[i+1] = table.concat(c, sep_ or ' ')
  end
  return table.concat(r, lsep_ or '\n')
end

function M.totable(x, r_)
  local nr, nc = x:sizes()
  local r = r_ or tbl_new(nr,0)
  assert(type(r) == 'table', "invalid argument, table expected")
  for i=0,nr-1 do
    local c = r[i+1] or tbl_new(nc,0)
    assert(type(c) == 'table', "invalid argument, table of tables expected")
    for j=0,nc-1 do c[j+1] = x:get0(i,j) end
    r[i+1] = c
  end
  return r
end

function M.fromtable (x, t)
  local nr, nc = x:sizes()
  assert(#t >= nr, "incompatible matrix-table column sizes")
  if type(t[1]) == 'table' then
    for i=0,nr-1 do
      local ti = t[i+1]
      assert(#ti >= nc, "incompatible matrix-table row sizes")
      for j=0,nc-1 do x:set0(i,j, ti[j+1]) end
    end
  else for i=0,nr-1 do x.data[i] = t[i+1] end
  end
  return x
end

M.__unm      = M.unm
M.__add      = M.add
M.__sub      = M.sub
M.__mul      = M.mul
M.__div      = M.div
M.__mod      = M.mod
M.__pow      = M.pow
M.__len      = M.size
M.__concat   = M.concat
M.__tostring = M.tostring

-- matrix-as-array behavior, unchecked bounds
M.__index    = function (self, idx)
  return isnum(idx) and self.data[idx] or M[idx]
end
M.__newindex = function (self, idx, val)
  self.data[idx] = val
end

ffi.metatype( 'matrix_t', M)
ffi.metatype('cmatrix_t', M)

------------------------------------------------------------------------------o
return matrix
