--[=[
 o----------------------------------------------------------------------------o
 |
 | Matrix constructor module
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
  - provides types and constructors to matrices

  Information:
  - this module is loaded by matrix and cmatrix, modules. It should
    not be loaded by users.

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  matctr -- linear algebra basics (types contructors)

SYNOPSIS
  This module should not be loaded directly, SEE ALSO.

DESCRIPTION
  The module matctr provides consistent definitions of matrices and
  complex matrices.

RETURN VALUES
  The constructors of matrices and complex matrices.

SEE ALSO
  math, gmath, complex, matrix, cmatrix
]]

-- modules -------------------------------------------------------------------o

local ffi   = require 'ffi'
local clib  = require 'cmad'
local gmath = require 'gmath'

-- ffi -----------------------------------------------------------------------o

ffi.cdef[[
typedef struct { int32_t nr, nc; double  data[?]; }  matrix_t;
typedef struct { int32_t nr, nc; complex data[?]; } cmatrix_t;
]]

-- threshold to use external allocator and save memory inside the 1GB limit
local mad_alloc = 1024

-- locals --------------------------------------------------------------------o

local istype  = ffi.istype

-- FFI type constructors
local  matrix_ctor = ffi.typeof( 'matrix_t')
local cmatrix_ctor = ffi.typeof('cmatrix_t')

-- implementation ------------------------------------------------------------o

function gmath.is_matrix (x)
  return type(x) == 'cdata' and istype('matrix_t', x)
end

function gmath.is_cmatrix (x)
  return type(x) == 'cdata' and istype('cmatrix_t', x)
end

function gmath.isa_matrix (x)
  return type(x) == 'cdata' and (istype('matrix_t', x) or istype('cmatrix_t', x))
end

local function matrix_alloc (nr, nc)
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

local function cmatrix_alloc (nr, nc)
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

local function fromtable (ctor, tbl)
  local nr, nc = #tbl or 1, type(tbl[1]) == 'table' and #tbl[1] or 1
  return ctor(nr,nc):fromtable(tbl)
end

local function matrix (nr, nc_)
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return fromtable(matrix_alloc, nr)
  elseif nr > 0 and nc > 0 then
    return matrix_alloc(nr, nc)
  else
    error("invalid argument to matrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

local function cmatrix (nr, nc_)
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return mat_fromtable(cmatrix_alloc, nr)
  elseif nr > 0 and nc > 0 then
    return cmatrix_alloc(nr, nc)
  else
    error("invalid argument to cmatrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

------------------------------------------------------------------------------o
return {
  matrix  = matrix,
  cmatrix = cmatrix,
}
