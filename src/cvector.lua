--[=[
 o----------------------------------------------------------------------------o
 |
 | Vector module (complex)
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
  - provides full set of functions and operations on complex vectors

  Information:
  - complex vectors are complex matrices [n x 1] implemented by module matrix

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  cvector

SYNOPSIS
  local cvector = require 'cvector'
  local v1 = cvector(3)                -- column cvector
  local v2 = cvector {1,2,3,4,5,6}     -- column cvector = {{1+0i},{2+0i},...}
  local v3 = v1:transpose()            -- row cvector, transpose conjugate

DESCRIPTION
  Complex vectors are complex matrices with one column.

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
  The constructor of complex vectors

SEE ALSO
  math, gmath, complex, vector, matrix, cmatrix
]]
 
-- modules -------------------------------------------------------------------o

local xmatrix = require 'xmatrix'

-- locals --------------------------------------------------------------------o

-- FFI type constructors
local cvector = xmatrix.cvector

-- implementation ------------------------------------------------------------o

-- implemented by the matrix module

------------------------------------------------------------------------------o
return cvector
