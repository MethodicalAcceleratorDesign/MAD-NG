--[=[
 o----------------------------------------------------------------------------o
 |
 | TPSA module (real)
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
  - real and complex TPSA are implemented by this module

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  tpsa

SYNOPSIS
  local tpsa = require 'tpsa'

DESCRIPTION
  The module tpsa implements the operators and math functions on TPSA:

RETURN VALUES
  The constructor of TPSA

SEE ALSO
  gmath, complex, matrix, cmatrix, ctpsa
]]
 
-- modules -------------------------------------------------------------------o

local ffi   = require 'ffi'
local clib  = require 'cmad'
local gmath = require 'gmath'
local xtpsa = require 'xtpsa'

local tbl_new = require 'table.new'

-- locals --------------------------------------------------------------------o

local isnum, iscpx, iscal, ismat, iscmat, isamat,
      real, imag, conj, ident, min,
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, mod, pow, tostring = 
      gmath.is_number, gmath.is_complex, gmath.is_scalar,
      gmath.is_matrix, gmath.is_cmatrix, gmath.isa_matrix,
      gmath.real, gmath.imag, gmath.conj, gmath.ident, gmath.min,
      gmath.abs, gmath.arg, gmath.exp, gmath.log, gmath.sqrt, gmath.proj,
      gmath.sin, gmath.cos, gmath.tan, gmath.sinh, gmath.cosh, gmath.tanh,
      gmath.asin, gmath.acos, gmath.atan, gmath.asinh, gmath.acosh, gmath.atanh,
      gmath.unm, gmath.mod, gmath.pow, gmath.tostring

local istype, cast, sizeof, fill = ffi.istype, ffi.cast, ffi.sizeof, ffi.fill

local cres = ffi.new 'complex[1]'

-- FFI type constructors

local tpsa  = xtpsa.tpsa
local ctpsa = xtpsa.ctpsa

-- implementation ------------------------------------------------------------o

ffi.metatype( 'tpsa_t', M)
ffi.metatype('ctpsa_t', M)

------------------------------------------------------------------------------o
return tpsa
