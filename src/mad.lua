--[=[
 o----------------------------------------------------------------------------o
 |
 | MAD environement (sandbox)
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
  - TODO
  
 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- module --------------------------------------------------------------------o

M.__help.self = [[
NAME
  mad -- Methodical Accelerator Design package

SYNOPSIS
  local mad = require "mad"

DESCRIPTION
  The MAD package provides all the modules and services required to run MAD.

RETURN VALUES
  The table of modules and services.

SEE ALSO
  None
]]

-- requires ------------------------------------------------------------------o

M.helper   = require "helper"
M.tester   = require "tester"

M.beam     = require "beam"
M.element  = require "element"
M.sequence = require "sequence"

-- end -----------------------------------------------------------------------o

return M
