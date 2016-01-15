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
typedef unsigned int bit_t; // mad_bit.h

struct desc { // warning: must be kept identical to C definition
  int   id;
  int   nmv, nv, nc;
  ord_t mo, ko, trunc;
};

struct tpsa { // warning: must be kept identical to C definition
  desc_t *desc;
  ord_t   lo, hi, mo;
  bit_t   nz;
  int     is_tmp;
  num_t   coef[?];
};

struct ctpsa { // warning: must be kept identical to C definition
  desc_t *desc;
  ord_t   lo, hi, mo;
  bit_t   nz;
  int     is_tmp;
  cnum_t  coef[?];
};
]]

-- locals --------------------------------------------------------------------o

local istype  = ffi.istype

-- FFI type constructors
local  tpsa_ctor   = ffi.typeof( 'tpsa_t')
local ctpsa_ctor   = ffi.typeof('ctpsa_t')

local strs_ctor    = ffi.typeof('str_t[?]')
local mono_ctor    = ffi.typeof('ord_t[?]')

-- threshold for external allocator
local mad_alloc = 1024

-- implementation ------------------------------------------------------------o

local function is_table(x)
  return type(x) == 'table'
end

local function is_desc(x)
  return type(x) == 'cdata' and istype('desc_t', x)
end

function gmath.is_tpsa (x)
  return type(x) == 'cdata' and istype('tpsa_t', x)
end

function gmath.is_ctpsa (x)
  return type(x) == 'cdata' and istype('ctpsa_t', x)
end

function gmath.isa_tpsa (x)
  return type(x) == 'cdata' and (istype('tpsa_t', x) or istype('ctpsa_t', x))
end

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

local isa_tpsa = gmath.isa_tpsa

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

-- nv: number of variables (if vo is a value)
-- vo: variables orders (array or value)
-- mo: map variables orders with mo[i] > vo[i]
-- nk: number of knobs (if ko is a value)
-- ko: knobs orders (array or value)
-- dk: max knobs 'cross' orders (degres)
-- ex0: {nv=2,vo=2 [, mo=3]}
-- ex1: {vo={2,2} [, mo={3,3}] [, v={'x', 'px'}] [, ko={1,1,1}] [, dk=2]}
-- ex2: {vo={2,2} [, mo={3,3}] [, v={'x', 'px'}] [, nk=3,ko=1] [, dk=2]}

local function desc (args)
  assert(args and args.vo, "not enough arguments for TPSA descriptor")

  local nv = args.nv or is_table(args.vo) and #args.vo or 0
  local nk = args.nk or is_table(args.ko) and #args.ko or 0
  local dk = args.dk or 0

  assert(nv > 0, "invalid number of variables")

  local cvar_ords =             mono_ctor(nv)
  local cmap_ords = args.mo and mono_ctor(nv) or nil
  local cvar_nams = args.v  and strs_ctor(nv) or nil
  local cknb_ords = nk > 0  and mono_ctor(nk) or nil

                    for i=1,nv do cvar_ords[i-1] = is_table(args.vo) and args.vo[i] or args.vo end
  if cmap_ords then for i=1,nv do cmap_ords[i-1] = is_table(args.mo) and args.mo[i] or args.mo end end
  if cvar_nams then for i=1,nv do cvar_nams[i-1] =                       args.v [i]            end end
  if cknb_ords then for i=1,nk do cknb_ords[i-1] = is_table(args.ko) and args.ko[i] or args.ko end end

  local desc

  if nk > 0 then
    desc = clib.mad_tpsa_desc_newk(nv, cvar_ords, cmap_ords, cvar_nams, nk, cknb_ords, dk)
  else
    desc = clib.mad_tpsa_desc_new (nv, cvar_ords, cmap_ords, cvar_nams)
  end

  tmp_stack[desc.id] = { top=0 }
  return desc
end

------------------------------------------------------------------------------o
return {
   tpsa =  tpsa,
  ctpsa = ctpsa,
   desc =  desc,
}
