--[=[
 o----------------------------------------------------------------------------o
 |
 | Test module
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
  tester -- run modules and functions tests

SYNOPSIS
  test = require "tester"
  test([mad_module_list])
  test.run(module)
  test.run(module.function)
  test.addModuleToTest(module_name)

DESCRIPTION
  The tester module runs the test of registered MAD modules and functions and
  return the statistics.
  
  test([mad_module_list])
    Runs the tests for all modules in the list. If those tests adds new
    dependencies, those are tested as well. An empty list will run all
    loaded modules.
  test.run(module/function)
    Runs the tests for a single module or function, plus dependencies.
  test.addModuleToTest(module_name)
    If a module depends on another module, run this function to add those
    dependencies to the list of modules to be tested.

RETURN VALUES
  The number of the tests failed and passed

SEE ALSO
  helper, module
]]

-- require -------------------------------------------------------------------o
-- local luaUnit = require "utest.luaUnit" -- to import...
local module  = require "module"

-- module --------------------------------------------------------------------o
local utest

function M.addModuleToTest( modname )
    utest:addModuleToTest(modname)
end

function M.run( mod, fun )
    if type(mod) == "table" then
        mod = module.get_module_name(mod)
    end
    if fun and type(fun) == "function" then
        fun = module.get_function_name(fun)
    end
    
    if type(mod) ~= "string" then
        error("Module isn't of type string")
    elseif fun and type(fun) ~= "string" then
        error("Function isn't of type string")
    end
    
    utest = luaUnit()
    if fun then
        utest:addFunctionToTest( mod, fun )
    else
        utest:addModuleToTest( mod )
    end
    
    utest:run()
end

-- metamethods ---------------------------------------------------------------o
local mt = {}; setmetatable(M, mt)

mt.__call = function (_, module_list)
    utest = luaUnit()
    for _, v in pairs(module_list) do
    	utest:addModuleToTest(v)
    end
    if #module_list == 0 then
        for k,v in pairs(module.get_all()) do
            utest:addModuleToTest(v)
        end
    end
    utest:run()
end

-- end -----------------------------------------------------------------------o

return M
