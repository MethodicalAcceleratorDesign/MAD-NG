--[=[
 o----------------------------------------------------------------------------o
 |
 | TPSA module (real)
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
  - provides full set of functions and operations on real TPSA

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  tpsa

SYNOPSIS
  local tpsa = require 'tpsa'

DESCRIPTION
  The module tpsa implements the operators and math functions on TPSA:

RETURN VALUES
  The constructor of TPSA

SEE ALSO
  gmath, complex, matrix, cmatrix, ctpsa
]]
 
-- modules -------------------------------------------------------------------o

local ffi   = require 'ffi'
local clib  = require 'cmad'
local gmath = require 'gmath'
local xtpsa = require 'xtpsa'

local tbl_new = require 'table.new'

-- locals --------------------------------------------------------------------o

local istype, cast, sizeof, fill = ffi.istype, ffi.cast, ffi.sizeof, ffi.fill
local min, max = math.min, math.max
local isnum, istable, istpsa = gmath.is_number, gmath.is_table, gmath.is_tpsa

-- FFI type constructors
local tpsa = xtpsa.tpsa
local desc = xtpsa.desc
local mono = xtpsa.mono

local cint_arr = ffi.typeof('const int[?]')

-- implementation ------------------------------------------------------------o


function M.mono (t, tbl_)
  return mono(tbl_ or t.d.nv)
end

-- initialization ------------------------------------------------------------o

function M.copy (t, r_)
  local r = r_ or t:tpsa()
  clib.mad_tpsa_copy(t, r)
  return r
end

M.clear  = clib.mad_tpsa_clear
M.scalar = clib.mad_tpsa_scalar

-- indexing ------------------------------------------------------------------o

function M.get_idx (t,tbl)
  local m = istable(tbl) and mono(tbl) or tbl
  return clib.mad_tpsa_midx(t, m.n, m.ord)
end

function M.get_idx_sp (t,tbl)
  -- tbl = {idx1, ord1, idx2, ord2, ... }
  local n = #tbl
  local m = cint_arr(n)
  for i=1,n do m[i-1] = tbl[i] end
  return clib.mad_tpsa_midx_sp(t, n, m)
end

function M.get_mono (t, i, r_)
  local  m = r_ or t:mono()
  return m, clib.mad_tpsa_mono(t, m.n, m.ord, i)
end

-- peek & poke ---------------------------------------------------------------o

M.get0 = clib.mad_tpsa_get0

function M.get (t, m)
  if isnum(m) then
    return clib.mad_tpsa_geti(t, m)
  end
  m = istable(m) and mono(m) or m
  return clib.mad_tpsa_getm(t, m.n, m.ord)
end

function M.get_sp (t, tbl)
  -- tbl = {idx1, ord1, idx2, ord2, ... }
  local n = #tbl
  local m = cint_arr(n)
  for i=1,n do m[i-1] = tbl[i] end
  return clib.mad_tpsa_getm_sp(t, n, m)
end

function M.set0 (t, a, b)
  if b == nil then a, b = 0, a end
  clib.mad_tpsa_set0(t,a,b)
end

function M.set (t, i, a, b)
  if b == nil then a, b = 0, a end
  if isnum(i) then
    clib.mad_tpsa_seti(t, i, a, b)
  end
  local m = istable(i) and mono(i) or i
  clib.mad_tpsa_setm(t, m.n, m.ord, a, b)
end

function M.set_sp (t, tbl, a, b)
  if b == nil then a, b = 0, a end
  -- tbl = {idx1, ord1, idx2, ord2, ... }
  local n = #tbl
  local m = cint_arr(n)
  for i=1,n do m[i-1] = tbl[i] end
  clib.mad_tpsa_setm_sp(t, n, m, a, b)
end

-- unary operators -----------------------------------------------------------o

function M.abs (t, r_)
  local r = r_ or t:tpsa()
  clib.mad_tpsa_abs(t,r)
end

M.nrm1 = clib.mad_tpsa_nrm1
M.nrm2 = clib.mad_tpsa_nrm2

function M.der (t, ivar, r_)
  local r = r_ or t:tpsa()
  clib.mad_tpsa_der(t, r, ivar)
  return r
end

function M.mder (t, tbl, r_)
  local r = r_ or t:tpsa()
  local m = mono(tbl)
  clib.mad_tpsa_mder(t, r, m.n, m.ord)
  return r
end

function M.scale (t, val, r_)
  local r = r_ or t:tpsa()
  clib.mad_tpsa_scl(t, val, r)
  return r
end

-- binary operators ----------------------------------------------------------o

function M.add (a, b, r_)
  local r
  if isnum(a) then       -- num + tpsa
    if b.hi == 0 then return a + b.coef[0] end
    r = b:copy(r_)
    clib.mad_tpsa_set0(r, 1, a)
  elseif isnum(b) then   -- tpsa + num
    if a.hi == 0 then return a.coef[0] + b end
    r = a:copy(r_)
    clib.mad_tpsa_set0(r, 1, b)
  elseif istpsa(b) then  -- tpsa + tpsa
    r = r_ or a:tpsa(max(a.mo,b.mo))
    clib.mad_tpsa_add(a, b, r)
  else error("invalid GTPSA (+) operands") end
  return r
end

function M.sub (a, b, r_)
  local r
  if isnum(a) then       -- num - tpsa
    if b.hi == 0 then return a - b.coef[0] end
    r = r_ or b:tpsa()
    clib.mad_tpsa_scl (b, -1, r)
    clib.mad_tpsa_set0(r,  1, a)
  elseif isnum(b) then   -- tpsa - num
    if a.hi == 0 then return a.coef[0] - b end
    r = a:copy(r_)
    clib.mad_tpsa_set0(r, 1, -b)
  elseif istpsa(b) then  -- tpsa - tpsa
    r = r_ or a:tpsa(max(a.mo,b.mo))
    clib.mad_tpsa_sub(a, b, r)
  else error("invalid GTPSA (-) operands") end
  return r
end

function M.mul (a, b, r_)
  local r
  if isnum(a) then       -- num * tpsa
    if b.hi == 0 then return a * b.coef[0] end
    r = r_ or b:tpsa()
    clib.mad_tpsa_scl(b, a, r)
  elseif isnum(b) then   -- tpsa * num
    if a.hi == 0 then return a.coef[0] * b end
    r = r_ or a:tpsa()
    clib.mad_tpsa_scl(a, b, r)
  elseif istpsa(b) then  -- tpsa * tpsa
    r = r_ or a:tpsa(max(a.mo,b.mo))
    clib.mad_tpsa_mul(a, b, r)
  else error("invalid GTPSA (*) operands") end
  return r
end

function M.div (a, b, r_)
  local r
  if isnum(a) then       -- num / tpsa
    if b.hi == 0 then return a / b.coef[0] end
    r = r_ or b:tpsa()
    clib.mad_tpsa_inv(b, a, r)
  elseif isnum(b) then   -- tpsa / num
    if a.hi == 0 then return a.coef[0] / b end
    r = r_ or a:tpsa()
    clib.mad_tpsa_scl(a, 1/b, r)
  elseif istpsa(b) then  -- tpsa / tpsa
    r = r_ or a:tpsa(max(a.mo,b.mo))
    clib.mad_tpsa_div(a, b, r)
  else error("invalid GTPSA (/) operands") end
  return r
end

function M.poisson(a, b, n, r_)
  local r = r_ or a:tpsa(max(a.mo,b.mo))
  clib.mad_tpsa_poisson(a, b, r, n)
  return r
end

------------------------------------------------------------------------------o

ffi.metatype('tpsa_t', M)

------------------------------------------------------------------------------o
return tpsa
