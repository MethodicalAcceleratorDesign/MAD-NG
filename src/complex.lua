--[=[
 o----------------------------------------------------------------------------o
 |
 | Complex number module
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
  - provides full set of functions and operations on complex numbers

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

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
    add, sub, mul, div,
    real, imag, conj,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    tostring.

RETURN VALUE
  The constructor of complex numbers

SEE ALSO
  math, gmath, cmatrix
]]

-- modules -------------------------------------------------------------------o

local ffi     = require 'ffi'
local clib    = require 'cmad'
local gmath   = require 'gmath'
local cmatrix = require 'cmatrix'

-- locals --------------------------------------------------------------------o

local isnum, iscpx, iscal, ismat, iscmat, tostring =
      gmath.is_number, gmath.is_complex, gmath.is_scalar,
      gmath.is_matrix, gmath.is_cmatrix, gmath.tostring

local istype = ffi.istype
local abs = math.abs

local cres = ffi.new 'complex[1]'

-- FFI type constructors
local complex = ffi.typeof 'complex'

-- implementation ------------------------------------------------------------o

function M.real  (x) return x.re end
function M.imag  (x) return x.im end
function M.conj  (x) return complex( x.re, -x.im) end

function M.abs   (x) return clib.mad_cnum_abs (x.re, x.im) end
function M.arg   (x) return clib.mad_cnum_arg (x.re, x.im) end

function M.exp   (x) clib.mad_cnum_exp   (x.re, x.im, cres) ; return cres[0] end
function M.log   (x) clib.mad_cnum_log   (x.re, x.im, cres) ; return cres[0] end
function M.sqrt  (x) clib.mad_cnum_sqrt  (x.re, x.im, cres) ; return cres[0] end
function M.proj  (x) clib.mad_cnum_proj  (x.re, x.im, cres) ; return cres[0] end

function M.sin   (x) clib.mad_cnum_sin   (x.re, x.im, cres) ; return cres[0] end
function M.cos   (x) clib.mad_cnum_cos   (x.re, x.im, cres) ; return cres[0] end
function M.tan   (x) clib.mad_cnum_tan   (x.re, x.im, cres) ; return cres[0] end
function M.sinh  (x) clib.mad_cnum_sinh  (x.re, x.im, cres) ; return cres[0] end
function M.cosh  (x) clib.mad_cnum_cosh  (x.re, x.im, cres) ; return cres[0] end
function M.tanh  (x) clib.mad_cnum_tanh  (x.re, x.im, cres) ; return cres[0] end

function M.asin  (x) clib.mad_cnum_asin  (x.re, x.im, cres) ; return cres[0] end
function M.acos  (x) clib.mad_cnum_acos  (x.re, x.im, cres) ; return cres[0] end
function M.atan  (x) clib.mad_cnum_atan  (x.re, x.im, cres) ; return cres[0] end
function M.asinh (x) clib.mad_cnum_asinh (x.re, x.im, cres) ; return cres[0] end
function M.acosh (x) clib.mad_cnum_acosh (x.re, x.im, cres) ; return cres[0] end
function M.atanh (x) clib.mad_cnum_atanh (x.re, x.im, cres) ; return cres[0] end

function M.__eq  (x, y)
  if iscal(y) then
    x, y = complex(x), complex(y)
    return x.re == y.re and x.im == y.im
  else
    return y == x
  end
end

function M.unm (x)
  return complex(-x.re, -x.im)
end

function M.add (x, y, r_)
  if isnum(x) then -- num + cpx
    return complex(x + y.re, y.im)
  elseif isnum(y) then -- cpx + num
    return complex(x.re + y, x.im)
  elseif iscpx(y) then -- cpx + cpx
    return complex(x.re + y.re, x.im + y.im)
  end

  -- iscpx(x)
  local r = r_ or cmatrix(y:sizes())
  assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")

  if ismat(y) then -- cpx + mat => vec + cpx
    clib.mad_vec_addc(y.data, x.re, x.im, r.data, r:size())
  elseif iscmat(y) then -- cpx + cmat => cvec + cpx
    clib.mad_cvec_addc(y.data, x.re, x.im, r.data, r:size())
  else error("incompatible complex (+) operands") end
  return r
end

function M.sub (x, y, r_)
  if isnum(x) then -- num - cpx
    return complex(x - y.re, - y.im)
  elseif isnum(y) then -- cpx - num
    return complex(x.re - y, x.im)
  elseif iscpx(y) then -- cpx - cpx
    return complex(x.re - y.re, x.im - y.im)
  end

  -- iscpx(x)
  local r = r_ or cmatrix(y:sizes())
  assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")

  if ismat(y) then -- cpx - mat => cpx - vec
    clib.mad_vec_subc(y.data, x.re, x.im, r.data, r:size())
  elseif iscmat(y) then -- cpx - cmat => cpx - cvec
    clib.mad_cvec_subc(y.data, x.re, x.im, r.data, r:size())
  else error("incompatible complex (-) operands") end
  return r
end

function M.mul (x, y, r_)
  if isnum(x) then -- num * cpx
    return complex(x * y.re, x * y.im)
  elseif isnum(y) then -- cpx * num
    return complex(x.re * y, x.im * y)
  elseif iscpx(y) then -- cpx * cpx
    return complex(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re)
  end

  -- iscpx(x)
  local r = r_ or cmatrix(y:sizes())
  assert(y:rows() == r:rows() and y:cols() == r:cols(), "incompatible cmatrix sizes")
    
  if ismat(y) then -- cpx * mat => vec * cpx
    clib.mad_vec_mulc(y.data, x.re, x.im, r.data, r:size())
  elseif iscmat(y) then -- cpx * cmat => cvec * cpx
    clib.mad_cvec_mulc(y.data, x.re, x.im, r.data, r:size())
  else error("invalid complex (*) operands") end
  return r
end

function M.div (x, y, r_, rcond_)
  if isnum(x) then -- num / cpx
    clib.mad_cnum_div(x, 0, y.re, y.im, cres)
    return cres[0] 
  elseif isnum(y) then -- cpx / num
    return complex(x.re / y, x.im / y)
  elseif iscpx(y) then -- cpx / cpx
    clib.mad_cnum_div(x.re, x.im, y.re, y.im, cres)
    return cres[0]
  end

  -- iscpx(x)
  local r = r_ or cmatrix(y:tsizes())
  assert(y:rows() == r:cols() and y:cols() == r:rows(), "incompatible cmatrix sizes")

  if ismat(y) then -- cpx / mat
    clib.mad_mat_invc (y.data, x.re, x.im, r.data, y:rows(), y:cols(), rcond_ or -1)
  elseif iscmat(y) then -- cpx / cmat
    clib.mad_cmat_invc(y.data, x.re, x.im, r.data, y:rows(), y:cols(), rcond_ or -1)
  else error("invalid complex (/) operands") end
  return r
end

function M.mod (x, y)
  error("cpx % cpx: NYI")
end

function M.pow (x, y)
  if isnum(y) then
    if y <  0 then x, y = 1/x, -y end
    if y == 2 then return x*x end
    if y == 1 then return x end
    if y == 0 then return 1 end
    y = complex(y)
  end

  if iscpx(y) then
    x = complex(x)
    clib.mad_cnum_pow(x.re, x.im, y.re, y.im, cres)
    return cres[0]
  end

  error("incompatible complex (^) operands")
end

function M.tostring (x)
      if x.im == 0 then return                        tostring(x.re)
  elseif x.re == 0 then return string.format('%si',                  tostring(x.im))
  elseif x.im <  0 then return string.format('%s%si', tostring(x.re),tostring(x.im))
  else                  return string.format('%s+%si',tostring(x.re),tostring(x.im))
  end
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

ffi.metatype('complex', M)

------------------------------------------------------------------------------o
return complex
