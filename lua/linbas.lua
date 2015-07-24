local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  linbas -- linear algebra basics (types contructors)

SYNOPSIS
  This module should not be loaded directly, SEE ALSO.

DESCRIPTION
  The module linbas provides consistent definitions of matrices and
  complex matrices.

RETURN VALUES
  The constructors of matrices and complex matrices.

SEE ALSO
  math, gmath, complex, matrix, cmatrix, linalg
]]

-- DEFS ------------------------------------------------------------------------

local ffi   = require 'ffi'
local gmath = require 'gmath'

ffi.cdef[[
typedef struct { int32_t nr, nc; double  data[?]; }  matrix_t;
typedef struct { int32_t nr, nc; complex data[?]; } cmatrix_t;
]]

local  mat_ct = ffi.typeof  'matrix_t'
local cmat_ct = ffi.typeof 'cmatrix_t'
local istype  = ffi.istype

function gmath.is_matrix (x)
  return type(x) == 'cdata' and istype('matrix_t', x)
end

function gmath.is_cmatrix (x)
  return type(x) == 'cdata' and istype('cmatrix_t', x)
end

function gmath.isa_matrix (x)
  return type(x) == 'cdata' and (istype('matrix_t', x) or istype('cmatrix_t', x))
end

local function mat_alloc (ct, nr, nc)
  local r = ct(nr*nc) -- default init: compiled for size <= 128 bytes
  r.nr, r.nc = nr, nc
  return r
end

local function mat_fromtable (ct, t)
  local nr, nc = #t, type(t[1]) == 'table' and #t[1] or 1
  local r = mat_alloc(ct, nr, nc)
  return r:fromtable(t)
end

local function matrix (nr, nc_)
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return mat_fromtable(mat_ct, nr)
  elseif nr > 0 and nc > 0 then
    return mat_alloc(mat_ct, nr, nc)
  else
    error("invalid argument to matrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

local function cmatrix (nr, nc_)
  local nc = nc_ or 1
  if type(nr) == 'table' then
    return mat_fromtable(cmat_ct, nr)
  elseif nr > 0 and nc > 0 then
    return mat_alloc(cmat_ct, nr, nc)
  else
    error("invalid argument to cmatrix constructor, expecting (rows[,cols]) or table [of tables]")
  end
end

-- END -------------------------------------------------------------------------
return {
  matrix  = matrix,
  cmatrix = cmatrix,
}
