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
    (minus) -, +, -, *, /, %, ^, ==, #, ..,
    unm, add, sub, mul, div, mod, pow,
    size, sizes, get, set, get0, set0,
    zeros, ones, unit, fill, copy,
    real, imag, conj, norm, angle,
    dot, inner, cross, mixed, outer, (products)
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, foreach, map, map2,
    concat, tostring, totable, fromtable.

REMARK
  check_bounds  =true checks out of bounds indexes in get , set
  check_bounds0 =true checks out of bounds indexes in get0, set0

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
