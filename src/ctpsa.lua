--[=[
 o----------------------------------------------------------------------------o
 |
 | TPSA module (complex)
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
  - provides full set of functions and operations on real and complex TPSA

  Information:
  - real and complex TPSA are implemented the module tpsa (not this one)

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  ctpsa

SYNOPSIS
  local ctpsa = require 'ctpsa'

DESCRIPTION
  The module ctpsa implements the operators and math functions on
  complex TPSA:
  TODO
    (minus) -, +, -, *, /, %, ^, ==, #, [], ..,
    unm, add, sub, mul, div, mod, pow, emul, ediv,
    rows, cols, size, sizes, get, set, get0, set0,
    zeros, ones, unit, fill, copy,
    get_row, get_col, get_diag, get_sub,
    set_row, set_col, set_diag, set_sub,
    transpose, t, trans, ctrans,
    real, imag, conj, norm, angle, trace, tr,
    dot, inner, cross, mixed, outer,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, foreach, map, map2, maps,
    concat, reshape, tostring, totable, fromtable,
    check_bounds.

REMARK:
  By default, check_bounds is true.

RETURN VALUES
  The constructor of complex TPSA

SEE ALSO
  gmath, complex, matrix, cmatrix, tpsa
]]
 
-- modules -------------------------------------------------------------------o

local xtpsa = require 'xtpsa'

-- locals --------------------------------------------------------------------o

-- FFI type constructors
local ctpsa = xtpsa.ctpsa

-- implementation ------------------------------------------------------------o

-- implemented by the matrix module

------------------------------------------------------------------------------o
return ctpsa
