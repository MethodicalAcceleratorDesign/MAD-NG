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
  - provides full set of functions and operations on real TPSA

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

local istype, cast, sizeof, fill = ffi.istype, ffi.cast, ffi.sizeof, ffi.fill

-- FFI type constructors
local desc = xtpsa.desc
local tpsa = xtpsa.tpsa

-- implementation ------------------------------------------------------------o

ffi.metatype( 'tpsa_t', M)

------------------------------------------------------------------------------o
return tpsa
