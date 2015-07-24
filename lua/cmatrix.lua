local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  cmatrix

SYNOPSIS
  local cmatrix = require 'cmatrix'
  local m1 = cmatrix(3)                 -- column cmatrix = cmatrix(3,1)
  local m2 = cmatrix(2,3)
  local m3 = cmatrix {{1,2+2i},{3,4+2i},{5,6+2i}}
  local m4 = cmatrix {1,2,3,4,5,6}      -- column cmatrix = {{1+0i},{2+0i},...}
  local m5 = cmatrix {{1,2,3,4,5,6}}    -- row cmatrix
  local m6 = m1:transpose()             -- row cmatrix, transpose conjugate

DESCRIPTION
  The module cmatrix implements the operators and math functions on
  complex matrices:
    (minus) -, +, -, *, /, %, ^, ==, #, [], ..,
    unm, add, sub, mul, div, mod, pow, emul, ediv,
    rows, cols, size, sizes, get, set, get0, set0,
    zeros, ones, unit, fill, copy, transpose, trans, t,
    get_row, get_col, get_diag, set_row, set_col, set_diag,
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

RELATIONS (3D GEOMETRY)
  inner prod:   u'.v = |u|.|v| cos(u^v)
  cross prod:   uxv = |u|.|v| sin(u^v) \vec{n}
  mixed prod:   (uxv)'.w = u'.(vxw) = det(u,v,w)
  outer prod:   u.v' = matrix
  dble xprod:   ux(vxw) = (u.w) \vec{v} - (u.v) \vec{w}
                (uxv)xw = (u.w) \vec{v} - (v.w) \vec{u}
  norm      :   |u| = sqrt(u'.u)
  angle     :   u^v = acos(u'.v / |u|.|v|)  in [0,pi] (or [-pi,pi] if n)
  unit      :   u / |u|
  projection:   u'.v
  projector :   I -   u.u' / u'.u
  reflector :   I - 2 u.u' / u'.u
  area      :   |uxv|
  volume    :   |(uxv)'.w|
  unitary   :   |u| = 1
  orthogonal:   u'.v = 0
  collinear :   |uxv| = 0
  coplanar  :   (uxv)'.w = 0

RETURN VALUES
  The constructor of complex matrices

SEE ALSO
  math, gmath, complex, matrix, cmatrix, linalg
]]
 
-- DEFS ------------------------------------------------------------------------

local linbas = require 'linbas'
local matrix = require 'matrix'

-- implemented by the matrix module

-- END -------------------------------------------------------------------------
return linbas.cmatrix
