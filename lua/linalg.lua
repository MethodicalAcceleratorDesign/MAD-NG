local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  linalg

SYNOPSIS
  local linalg = require 'linalg'
  local vector, cvector = linalg.vector, linalg.cvector
  local matrix, cmatrix = linalg.matrix, linalg.cmatrix

DESCRIPTION
  The module linalg provides consistent definitions of vectors, complex vectors,
  matrices and complex matrices.

RETURN VALUES
  The constructors of vectors, complex vectors, matrices and complex matrices.

SEE ALSO
  math, gmath, complex, vector, cvector, matrix, cmatrix
]]

-- DEFS ------------------------------------------------------------------------

local ffi  = require 'ffi'
local cmad = require 'cmad'
local gm   = require 'gmath'

ffi.cdef[[
typedef struct { int32_t n;      double  data[?]; }  vector_t;
typedef struct { int32_t n;      complex data[?]; } cvector_t;
typedef struct { int32_t nr, nc; double  data[?]; }  matrix_t;
typedef struct { int32_t nr, nc; complex data[?]; } cmatrix_t;
]]

-- locals
local  vec_ct = ffi.typeof  'vector_t'
local cvec_ct = ffi.typeof 'cvector_t'
local  mat_ct = ffi.typeof  'matrix_t'
local cmat_ct = ffi.typeof 'cmatrix_t'
local istype  = ffi.istype

-- extend gmath
function gm.is_vector  (x) return type(x) == 'cdata' and istype( 'vector_t', x) end
function gm.is_cvector (x) return type(x) == 'cdata' and istype('cvector_t', x) end
function gm.is_matrix  (x) return type(x) == 'cdata' and istype( 'matrix_t', x) end
function gm.is_cmatrix (x) return type(x) == 'cdata' and istype('cmatrix_t', x) end

-- vec --------------------------------

local function vec_alloc (ct, n)
  local r = ct(n) -- default init: compiled for size <= 128 bytes
  r.n = n
  return r
end

local function vector_from_table(ct, tbl)
  local r = vec_alloc(ct, #tbl)
  return r:set_table(tbl)
end

local function vector (n)
  if type(n) == 'table' then
    return vector_from_table(vec_ct, n)
  elseif n > 0 then
    return vec_alloc(vec_ct, n)
  else
    error("invalid argument to vector constructor, expecting size or table")
  end
end

local function cvector (n)
  if type(n) == 'table' then
    return vector_from_table(cvec_ct, n)
  elseif n > 0 then
    return vec_alloc(cvec_ct, n)
  else
    error("invalid argument to cvector constructor, expecting size or table")
  end
end

-- mat --------------------------------

local function mat_alloc (ct, nr, nc)
  local r = ct(nr*nc) -- default init: compiled for size <= 128 bytes
  r.nr, r.nc = nr, nc
  return r
end

local function matrix_from_table (ct, tbl)
  local r = mat_alloc(ct, #tbl, #tbl[1])
  return r:set_table(tbl)
end

local function matrix (nr, nc)
  if type(nr) == 'table' then
    return matrix_from_table(mat_ct, nr)
  elseif nr > 0 and nc > 0 then
    return mat_alloc(mat_ct, nr, nc)
  else
    error("invalid argument to matrix constructor, expecting (rows,cols) or table of tables")
  end
end

local function cmatrix (nr, nc)
  if type(nr) == 'table' then
    return matrix_from_table(cmat_ct, nr)
  elseif nr > 0 and nc > 0 then
    return mat_alloc(cmat_ct, nr, nc)
  else
    error("invalid argument to cmatrix constructor, expecting (rows,cols) or table of tables")
  end
end

-- END -------------------------------------------------------------------------
return {
  cmad    = cmad,
  vector  = vector,
  matrix  = matrix,
  cvector = cvector,
  cmatrix = cmatrix,
}
