--[=[
 o----------------------------------------------------------------------------o
 |
 | Module module (for helper)
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

M.__help.self = [=[
NAME
  module -- register MAD modules and functions

SYNOPSIS
  module = require "module"
  mod_name           = module.get_module_name(mod)
  fun_name, mod_name = module.get_function_name(fun)
  mod_list, fun_list = module.get_all()

DESCRIPTION
  The module module register MAD modules and functions for name lookup. MAD
  modules must have '__help' and '__test' tables.

  A MAD module should typically start by:

    local M = { __help={}, __test={} }
    M.__help.self = [[
      module documentation
    ]]
    M.__help.func_name = [[
      function documentation
    ]]

RETURN VALUES
  The interface to the MAD modules database.

ERRORS
  If the module does not contain a '__help' and '__test' tables, an invalid
  argument error is raised.

SEE ALSO
  helper, tester, object
]=]

-- locals ----------------------------------------------------------------------

local registered = false
local registered_module = {}    -- { [mod] = 'mod_name' }
local registered_function = {}  -- { [fun] = { fun_name = 'fun_name', mod_name = 'mod_name' } }

-- functions -------------------------------------------------------------------

local register_function = function (mod, mod_name)
  for fun_name,fun in pairs(mod) do
    if type(fun) == 'function' then
      registered_function[fun] = { fun_name = fun_name, mod_name = mod_name }
    end
  end
end

local register_module = function (mod_name)
  local mod = package.loaded[mod_name]

  if registered_module[mod] ~= nil then return end

  if not mod.__help then
    error(("module '%s' has NO help"):format(mod_name))
  end
  if not mod.__test then
    error(("module '%s' has NO test"):format(mod_name))
  end

  register_function(mod, mod_name)
  registered_module[mod] = mod_name
end

local register_mad_modules = function ()
  local mad = require 'mad'

  for mod_name,mod in pairs(package.loaded) do
    if type(mod) == 'table' and mod.__help and mod.__test then
      register_module(mod_name)
    end
  end

  registered = true
end

-- methods -------------------------------------------------------------------o

M.reset = function ()
  registered = false
  registered_module = {}
  registered_function = {}
end

M.get_all = function ()
  if not registered then
    register_mad_modules()
  end

  return registered_module, registered_function
end

M.get_module_name = function (mod)
  if not registered then
    register_mad_modules()
  end

  return registered_module[mod]
end

M.get_function_name = function (fun)
  if not registered then
    register_mad_modules()
  end
    
  local info = registered_function[fun]
  if info then
    return info.fun_name, info.mod_name
  else
    return nil, nil
  end
end

-- end -----------------------------------------------------------------------o

return M
