local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  cvector

SYNOPSIS
  local cvector = require 'cvector'
  local v1 = cvector {1,2+2i,3,4+2i,5,6+2i}
  local v2 = cvector(6)
  local r1 = v1:dot(v2)
  local v3 = v1:cross(v2)

DESCRIPTION
  The module cvector implements the operators and math functions on
  complex vectors:
    (minus) -, +, -, *, /, %, ^, ==,
    unm, add, sub, mul, div, mod, pow,
    size, sizes, get, set, get0, set0,
    zeros, ones, fill, copy, set_table,
    real, imag, conj, norm, angle,
    inner_prod, outer_prod, cross_prod,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, map, map2, tostring, totable.

REMARK
  check_bounds  can be set to true to check out of bounds indexes in get , set
  check_bounds0 can be set to true to check out of bounds indexes in get0, set0

RETURN VALUES
  The constructor of complex vectors

SEE ALSO
  math, gmath, complex, vector, matrix, cmatrix
]]
 
-- DEFS ------------------------------------------------------------------------

local linalg = require 'linalg'
local vector = require 'vector'

-- implemented by the vector module

-- END -------------------------------------------------------------------------
return linalg.cvector
