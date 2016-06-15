--[=[
 o-----------------------------------------------------------------------------o
 |
 | Lua core feature regression tests
 |
 | Methodical Accelerator Design - Copyright CERN 2015+
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
  - Provide a small set of test suites for some Lua features.

 o-----------------------------------------------------------------------------o
]=]

-- locals ---------------------------------------------------------------------o

local luaunit = require 'luaunit'
local assertEquals, assertAlmostEquals =
      luaunit.assertEquals, luaunit.assertAlmostEquals

-- regression test suite ------------------------------------------------------o

TestLuaCore = {}

local function isPrimeDivisible(self, c)
  for i=3, self.prime_count do
    if self.primes[i] * self.primes[i] > c then break end
    if c % self.primes[i] == 0 then return true end
  end
  return false
end

local function addPrime(self, c)
  self.prime_count = self.prime_count + 1
  self.primes[self.prime_count] = c
end

local function getPrimes(n)
  local p = { prime_count=3, primes={ 1,2,3 } }
  local c = 5
  while p.prime_count < n do
    if not isPrimeDivisible(p, c) then
      addPrime(p, c)
    end
    c = c + 2 -- process only odd numbers
  end
  return p
end

function TestLuaCore:testPrimes()
  local p = getPrimes(1e3)
  assertEquals( p.primes[p.prime_count-0], 7907 )
  assertEquals( p.primes[p.prime_count-1], 7901 )
  assertEquals( p.primes[p.prime_count-2], 7883 )
  assertEquals( p.primes[p.prime_count-3], 7879 )
  assertEquals( p.primes[p.prime_count-4], 7877 )
  assertEquals( p.primes[p.prime_count-5], 7873 )
end

local function find_duplicates(inp)
  local res = {}
  for _,v in ipairs(inp) do
    res[v] = res[v] and res[v]+1 or 1
  end
  for _,v in ipairs(inp) do
    if res[v] and res[v] > 1 then
      res[#res+1] = v
    end
    res[v] = nil
  end
  return res
end

local function find_duplicates2(inp, res)
  for i=1,#res do res[i]=nil end
  for _,v in ipairs(inp) do
    res[v] = res[v] and res[v]+1 or 1
  end
  for _,v in ipairs(inp) do
    if res[v] and res[v] > 1 then
      res[#res+1] = v
    end
    res[v] = nil
  end
end

function TestLuaCore:testDuplicates()
  local inp = {'b','a','c','c','e','a','c','d','c','d'}
  local out = {'a','c','d'}
  local res = find_duplicates2(inp, res)
  assertEquals( find_duplicates(inp), out )
  assertEquals( res, out )
end

-- performance test suite -----------------------------------------------------o

Test_LuaCore = {}

function Test_LuaCore:testPrimes()
  local t0, p = os.clock()
  p = getPrimes(2e5)
  local dt = os.clock() - t0
  assertAlmostEquals( dt , 0.5, 1 )
  assertEquals( p.primes[p.prime_count], 2750131 )
end

function Test_LuaCore:testDuplicates()
  local inp = {'b','a','c','c','e','a','c','d','c','d'}
  local out = {'a','c','d'}
  local t0, res = os.clock()
  for i=1,5e5 do res = find_duplicates(inp) end
  local dt = os.clock() - t0
  assertAlmostEquals( dt , 0.4, 1 )
  assertEquals( res, out )
end

function Test_LuaCore:testDuplicates2()
  local inp = {'b','a','c','c','e','a','c','d','c','d'}
  local out = {'a','c','d'}
  local res = {}
  local t0 = os.clock()
  for i=1,5e5 do find_duplicates2(inp, res) end
  local dt = os.clock() - t0
  assertAlmostEquals( dt , 0.2, 1 )
  assertEquals( res, out )
end

-- end ------------------------------------------------------------------------o
