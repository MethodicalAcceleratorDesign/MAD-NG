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
    (minus) -, +, -, *, /, %, ^, ==,
    rows, cols, size, sizes, get, set, zeros, ones,
    unm, add, sub, mul, div, mod, pow, schur_prod,
    get_row, get_col, get_diag, transpose,
    set_row, set_col, set_diag, set_table,
    real, imag, conj, trace, norm, angle, inner_prod,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    copy, foldl, foldr, foreach, map, map2,
    tostring, totable.

RETURN VALUES
  The constructor of complex matrices

SEE ALSO
  math, gmath, complex, vector, cvector, matrix
]]
 
-- DEFS ------------------------------------------------------------------------

local linalg = require 'linalg'
local matrix = require 'matrix'

-- implemented in matrix module

-- END -------------------------------------------------------------------------
return linalg.cmatrix
