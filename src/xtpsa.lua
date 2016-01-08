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
typedef struct { // warning: must be kept identical to C definition
  int   id;
  int   nmv, nv, nc;
  ord_t mo, ko, trunc;
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

local function is_desc(x)
  return type(x) == 'cdata' and istype('desc_t', x)
end

local isa_tpsa = gmath.isa_tpsa

local function tpsa_alloc (desc, mo)
  local len, tpsa = desc.nc, nil
  if len < mad_alloc then
    tpsa = tpsa_ctor(len)
  else
    local siz = ffi.sizeof('tpsa_t', len)
    local ptr = clib.mad_malloc(siz)
    tpsa = ffi.gc(ffi.cast('tpsa_t&', ptr), clib.mad_free)
  end
  tpsa.desc = desc
  tpsa.lo, tpsa.hi, tpsa.mo = mo, 0, mo
  tpsa.nz, tpsa.is_tmp, tpsa.coef[0] = 0, 0, 0
  return tpsa
end

local function ctpsa_alloc (desc, mo)
  local len, tpsa = desc.nc, nil
  if len < (mad_alloc/2) then
    tpsa = tpsa_ctor(len)
  else
    local siz = ffi.sizeof('ctpsa_t', len)
    local ptr = clib.mad_malloc(siz)
    tpsa = ffi.gc(ffi.cast('ctpsa_t&', ptr), clib.mad_free)
  end
  tpsa.desc = desc
  tpsa.lo, tpsa.hi, tpsa.mo = mo, 0, mo
  tpsa.nz, tpsa.is_tmp, tpsa.coef[0] = 0, 0, 0
  return tpsa
end

local function tpsa (t, mo_)
  if isa_tpsa(t) then
    return tpsa_alloc(t.d, mo_ or t.mo)
  elseif is_desc(t) then
    return tpsa_alloc(  d, mo_ or d.mo)
  else
    error("invalid argument to tpsa constructor, expecting TPSA or descriptor")
  end
end

local function ctpsa (t, mo_)
  if isa_tpsa(t) then
    return ctpsa_alloc(t.d, mo_ or t.mo)
  elseif is_desc(t) then
    return ctpsa_alloc(  d, mo_ or d.mo)
  else
    error("invalid argument to ctpsa constructor, expecting TPSA or descriptor")
  end
end

------------------------------------------------------------------------------o
return {
   tpsa =  tpsa,
  ctpsa = ctpsa,
}
