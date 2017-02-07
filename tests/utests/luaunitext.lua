--[=[
 o-----------------------------------------------------------------------------o
 |
 | Luaunit extension for test suites
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - Provide some extensions to luaunit

 o-----------------------------------------------------------------------------o
]=]

local lu = require 'luaunit'
local assertEquals, assertAlmostEquals = lu.assertEquals, lu.assertAlmostEquals

lu.assertAllEquals = function (actual, expected, margin)
  assert(type(actual) == 'table', "invalid argument #1 (table expected)")

  if type(expected) ~= 'table' then
    local value = expected
    expected = {}
    for k in pairs(actual) do
      expected[k] = value
    end
  end

  if type(margin) ~= 'table' then
    local value = margin or 0
    margin = {}
    for k in pairs(actual) do
      margin[k] = value
    end
  end

  for k,v in pairs(actual) do
    if type(v) == 'number' and margin[k] ~= 0 then
      assertAlmostEquals (v, expected[k], margin[k])
    else
      assertEquals (v, expected[k])
    end
  end
end

-- tests ----------------------------------------------------------------------o

local assertAllEquals = lu.assertAllEquals
local assertErrorMsgContains = lu.assertErrorMsgContains

TestLuaUnitExt = {}

function TestLuaUnitExt:testAllEqualError()
  local errmsg = "invalid argument #1 (table expected)"
  local actual, expected = 0, 0
  assertErrorMsgContains(errmsg, assertAllEquals, actual, expected)
end

function TestLuaUnitExt:testAllEqualNumbers()
  local actual   = { -1.1*2, 0*2, 1.1*2 }
  local expected = { -2.2  , 0  , 2.2   }
  assertAllEquals (actual, expected)
end

function TestLuaUnitExt:testAllEqualNumbers2()
  local actual   = { -1.1-1e-15, -1.1, 0-1e-15, 0, 0+1e-15, 1.1, 1.1+1e-15 }
  local expected = { -1.1      , -1.1, 0      , 0, 0      , 1.1, 1.1       }
  local margin   = {      2e-15,    0,   1e-15, 0,   1e-15,   0,     2e-15 }
  assertAllEquals (actual, expected, margin)
end

function TestLuaUnitExt:testAllEqualMixed()
  local actual   = { -1.1-1e-15, {}, 0-1e-15, 0, 0+1e-15, {}, 1.1+1e-15 }
  local expected = { -1.1      , {}, 0      , 0, 0      , {}, 1.1       }
  local margin   = {      2e-15, {},   1e-15, 0,   1e-15, {},     2e-15 }
  assertAllEquals (actual, expected, margin)
end

function TestLuaUnitExt:testAllEqualKeyNumbers()
  local actual   = { x=-1.1*2, y=0*2, z=1.1*2 }
  local expected = { x=-2.2  , y=0  , z=2.2   }
  assertAllEquals (actual, expected)
end

function TestLuaUnitExt:testAllEqualKeyNumbers2()
  local actual   = { a=-1.1-1e-15, b=-1.1, c=0-1e-15, d=0, e=0+1e-15, f=1.1, g=1.1+1e-15 }
  local expected = { a=-1.1      , b=-1.1, c=0      , d=0, e=0      , f=1.1, g=1.1       }
  local margin   = { a=     2e-15, b=   0, c=  1e-15, d=0, e=  1e-15, f=  0, g=    2e-15 }
  assertAllEquals (actual, expected, margin)
end

function TestLuaUnitExt:testAllEqualKeyMixed()
  local actual   = { a=-1.1-1e-15, b={}, c=0-1e-15, d=0, e=0+1e-15, f={}, g=1.1+1e-15 }
  local expected = { a=-1.1      , b={}, c=0      , d=0, e=0      , f={}, g=1.1       }
  local margin   = { a=     2e-15, b={}, c=  1e-15, d=0, e=  1e-15, f={}, g=    2e-15 }
  assertAllEquals (actual, expected, margin)
end

function TestLuaUnitExt:testAllEqualMixedMixed()
  local actual   = { a=-1.1-1e-15, b={}, 0-1e-15, d=0, e=0+1e-15, {}, g=1.1+1e-15 }
  local expected = { a=-1.1      , b={}, 0      , d=0, e=0      , {}, g=1.1       }
  local margin   = { a=     2e-15, b={},   1e-15, d=0, e=  1e-15, {}, g=    2e-15 }
  assertAllEquals (actual, expected, margin)
end
