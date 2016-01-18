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

-- FFI type constructors
local desc = xtpsa.desc
local tpsa = xtpsa.tpsa

-- implementation ------------------------------------------------------------o

function M.copy(x, r_)
  local r = r_ or x:tpsa()
  clib.mad_tpsa_copy(x, r)
  return r
end

-- indexing -------------------------------------------------------------

function M.get_idx(t,m)
  return clib.mad_tpsa_midx(t,#m,mono_t(#m,m))
end

function M.get_idx_sp(t,m)
  return clib.mad_tpsa_midx_sp(t,#m,smono_t(#m,m))
end

function M.get_mono(t,i)
  local nv, ord = int_ptr(), ord_ptr()
  local cmono = clib.mad_tpsa_mono(t,i,nv,ord)
  local m = {}
  for i=1,nv[0] do
    m[i] = cmono[i-1]
  end
  return m, ord[0]
end

-- PEEK & POKE -----------------------------------------------------------------

M.clear  = clib.mad_tpsa_clear
M.get0   = clib.mad_tpsa_get0
M.get_at = clib.mad_tpsa_geti
M.scalar = clib.mad_tpsa_scalar

function M.get(t, m)
  return clib.mad_tpsa_getm(t, #m, mono_t(#m,m))
end

function M.get_sp(t,m)
  -- m = {idx1, ord1, idx2, ord2, ... }
  return clib.mad_tpsa_getm_sp(t, #m, smono_t(#m,m))
end

function M.set0(t, a,b)
  if not b then a, b = 0, a end
  clib.mad_tpsa_set0(t,a,b)
end


function M.set_at(t, i, a,b)
  if not b then a, b = 0, a end
  clib.mad_tpsa_seti(t,i,a,b)
end

function M.set(t, m, a,b)
  if not b then a, b = 0, a end
  clib.mad_tpsa_setm(t, #m, mono_t(#m, m), a, b)
end

function M.set_sp(t, m, a,b)
  if not b then a, b = 0, a end
  clib.mad_tpsa_setm_sp(t, #m, smono_t(#m,m), a, b)
end



ffi.metatype('tpsa_t', M)

------------------------------------------------------------------------------o
return tpsa
