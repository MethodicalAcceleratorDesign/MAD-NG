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

local ffi = require 'ffi'
local gm  = require 'gmath'

ffi.cdef[[
typedef struct { int32_t n;      double  data[?]; }  vector_t;
typedef struct { int32_t n;      complex data[?]; } cvector_t;
typedef struct { int32_t nr, nc; double  data[?]; }  matrix_t;
typedef struct { int32_t nr, nc; complex data[?]; } cmatrix_t;
]]

local  vec_ct = ffi.typeof  'vector_t'
local cvec_ct = ffi.typeof 'cvector_t'
local  mat_ct = ffi.typeof  'matrix_t'
local cmat_ct = ffi.typeof 'cmatrix_t'

-- extend gmath
local istype = ffi.istype
function gm.is_vector  (x) return type(x) == 'cdata' and istype( 'vector_t', x) end
function gm.is_cvector (x) return type(x) == 'cdata' and istype('cvector_t', x) end
function gm.is_matrix  (x) return type(x) == 'cdata' and istype( 'matrix_t', x) end
function gm.is_cmatrix (x) return type(x) == 'cdata' and istype('cmatrix_t', x) end

-- vec

local function vec_alloc (n)
  local r = vec_ct(n) -- default init: compiled for size <= 128 bytes
  r.n = n
  return r
end

local function vector_from_table(tbl)
  local r = vec_alloc(#tbl)
  return r:set_table(tbl)
end

local function vector (n)
  if type(n) == 'table' then
    return vector_from_table(n)
  end

  if n > 0 then
    return vec_alloc(n)
  end

  error("invalid argument to vector constructor, expecting size or table")
end

-- cvec

local function cvec_alloc (n)
  local r = cvec_ct(n) -- default init: compiled for size <= 128 bytes
  r.n = n
  return r
end

local function cvector_from_table (tbl)
  local r = cvec_alloc(#tbl)
  return r:set_table(tbl)
end

local function cvector (n)
  if type(n) == 'table' then
    return cvector_from_table(n)
  end

  if n > 0 then
    return cvec_alloc(n)
  end

  error("invalid argument to cvector constructor, expecting size or table")
end

-- mat

local function mat_alloc (nr, nc)
  local r = mat_ct(nr*nc) -- default init: compiled for size <= 128 bytes
  r.nr, r.nc = nr, nc
  return r
end

local function matrix_from_table (tbl)
  local r = mat_alloc(#tbl, #tbl[1])
  return r:set_table(tbl)
end

local function matrix (nr, nc)
  if type(nr) == 'table' then
    return matrix_from_table(nr)
  end

  if nr > 0 and nc > 0 then
    return mat_alloc(nr, nc)
  end

  error("invalid argument to matrix constructor, expecting (rows,cols) or table of tables")
end

-- cmat

local function cmat_alloc (nr, nc)
  local r = cmat_ct(nr*nc) -- default init: compiled for size <= 128 bytes
  r.nr, r.nc = nr, nc
  return r
end

local function cmatrix_from_table (tbl)
  local r = cmat_alloc(#tbl, #tbl[1])
  return r:set_table(tbl)
end

local function cmatrix (nr, nc)
  if type(nr) == 'table' then
    return cmatrix_from_table(nr)
  end

  if nr > 0 and nc > 0 then
    return cmat_alloc(nr, nc)
  end

  error("invalid argument to cmatrix constructor, expecting (rows,cols) or table of tables")
end

-- END -------------------------------------------------------------------------
return {
  vector  = vector,
  matrix  = matrix,
  cvector = cvector,
  cmatrix = cmatrix,
}


