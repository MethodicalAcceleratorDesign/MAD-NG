--[=[
 o----------------------------------------------------------------------------o
 |
 | TPSA constructor module
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
  - provides types and constructors to tpsa

  Information:
  - this module is loaded by tpsa and ctpsa modules. It should not be loaded
    by users.

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  xtpsa -- TPSA contructors

SYNOPSIS
  This module should not be loaded directly, SEE ALSO.

DESCRIPTION
  The module xtpsa provides consistent definitions of TPSA and complex TPSA.

RETURN VALUES
  The constructors of TPSA and complex TPSA.

SEE ALSO
  tpsa, ctpsa
]]

-- modules -------------------------------------------------------------------o

local ffi   = require 'ffi'
local clib  = require 'cmad'
local gmath = require 'gmath'

-- ffi -----------------------------------------------------------------------o

ffi.cdef[[
typedef struct {
  int id;     // warning: rest of fields are not exposed
} desc_t;

typedef struct { // warning: must be kept identical to C definition
  desc_t *desc;
  ord_t   lo, hi, mo;
  bit_t   nz;
  int     is_tmp;
  num_t   coef[?];
} tpsa_t;

typedef struct { // warning: must be kept identical to C definition
  desc_t *desc;
  ord_t   lo, hi, mo;
  bit_t   nz;
  int     is_tmp;
  cnum_t  coef[?];
} ctpsa_t;
]]

-- threshold to use external allocator and save memory inside the 1GB limit
local mad_alloc = 1024

-- locals --------------------------------------------------------------------o

local istype  = ffi.istype

-- FFI type constructors
local  tpsa_ctor = ffi.typeof( 'tpsa_t')
local ctpsa_ctor = ffi.typeof('ctpsa_t')

-- implementation ------------------------------------------------------------o

function gmath.is_tpsa (x)
  return type(x) == 'cdata' and istype('tpsa_t', x)
end

function gmath.is_ctpsa (x)
  return type(x) == 'cdata' and istype('ctpsa_t', x)
end

function gmath.isa_tpsa (x)
  return type(x) == 'cdata' and (istype('tpsa_t', x) or istype('ctpsa_t', x))
end

local function tpsa_alloc (nr, nc)
  -- TODO
  local len = nr*nc
  local mat
  if len < mad_alloc then
    mat = matrix_ctor(len)
  else
    local siz = ffi.sizeof('matrix_t', len)
    local ptr = clib.mad_malloc(siz)
    local mat = ffi.gc(ffi.cast('matrix_t&', ptr), clib.mad_free)
    ffi.fill(mat, siz)
  end
  mat.nr, mat.nc = nr, nc
  return mat
end

local function ctpsa_alloc (nr, nc)
  -- TODO
  local len = nr*nc
  local mat
  if len < (mad_alloc/2) then
    mat = cmatrix_ctor(len)
  else
    local siz = ffi.sizeof('cmatrix_t', len)
    local ptr = clib.mad_malloc(siz)
    local mat = ffi.gc(ffi.cast('cmatrix_t&', ptr), clib.mad_free)
    ffi.fill(mat, siz)
  end
  mat.nr, mat.nc = nr, nc
  return mat
end

local function tpsa (nr, nc_)
  -- TODO
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return fromtable(tpsa_alloc, nr)
  elseif nr > 0 and nc > 0 then
    return tpsa_alloc(nr, nc)
  else
    error("invalid argument to matrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

local function ctpsa (nr, nc_)
  -- TODO
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return mat_fromtable(cmatrix_alloc, nr)
  elseif nr > 0 and nc > 0 then
    return ctpsa_alloc(nr, nc)
  else
    error("invalid argument to cmatrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

------------------------------------------------------------------------------o
return {
   tpsa =  tpsa,
  ctpsa = ctpsa,
}
