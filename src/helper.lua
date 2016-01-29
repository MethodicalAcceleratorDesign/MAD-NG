--[=[
 o----------------------------------------------------------------------------o
 |
 | Help module
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

-- module ----------------------------------------------------------------------

M.__help.self = [[
NAME
  helper -- display modules and functions help on the console

SYNOPSIS
  help = require "helper"
  help(module)
  help(module.function)

DESCRIPTION
  The helper module displays the help of registered MAD modules and functions
  on the terminal.

RETURN VALUES
  True if the help was found, false otherwise.

SEE ALSO
  tester, module
]]

-- requires ------------------------------------------------------------------o

local module = require "module"

-- metamethods ---------------------------------------------------------------o

local MT = {}; setmetatable(M, MT)

function MT:__call(a)
  if type(a) == "table" then
    local mod = a
    local mod_name = module.get_module_name(mod)
    if mod_name and mod.__help.self then
      io.write(mod.__help.self)
      return true
    end

  elseif type(a) == "function" then
    local fun = a
    local fun_name, mod = module.get_function_name(fun)
    if fun_name and mod.__help[fun_name] then
      io.write(mod.__help[fun_name])
      return true
    end

  else
    io.write("No help found for ", a)
    return false
  end
end

-- end -----------------------------------------------------------------------o

return M
