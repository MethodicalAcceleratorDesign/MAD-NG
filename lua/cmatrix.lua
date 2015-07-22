local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  cmatrix

SYNOPSIS
  local cmatrix = require 'cmatrix'
  local m1 = cmatrix {{1,2+2i},{3,4+2i},{5,6+2i}}
  local m2 = cmatrix(2,3)
  local m3 = m1:transpose()
  m1:transpose(m2) 

DESCRIPTION
  The module cmatrix implements the operators and math functions on
  complex matrices:
    (minus) -, +, -, *, /, %, ^, ==, #, ..,
    unm, add, sub, mul, div, mod, pow, schur,
    rows, cols, size, sizes, get, set, get0, set0,
    zeros, ones, unit, fill, copy, transpose, t,
    get_row, get_col, get_diag, set_row, set_col, set_diag,
    real, imag, conj, norm, angle, dot, inner, trace, tr,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, foreach, map, map2,
    concat, reshape, tostring, totable, tovector, fromtable.

REMARK
  check_bounds  =true checks out of bounds indexes in get , set
  check_bounds0 =true checks out of bounds indexes in get0, set0

RETURN VALUES
  The constructor of complex matrices

SEE ALSO
  math, gmath, complex, vector, cvector, matrix
]]
 
-- DEFS ------------------------------------------------------------------------

local linalg = require 'linalg'
local matrix = require 'matrix'

-- implemented by the matrix module

-- END -------------------------------------------------------------------------
return linalg.cmatrix
