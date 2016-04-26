local Object = require 'object'

local Point = Object 'Point' {}

local p1 = Point      { x=1, y=2  }
local p2 = Point 'p2' { x=2, y=-1 }
local p3 = p2    'p3' { x=3  }

print(p1.name, p1.x, p1.y)
print(p2.name, p2.x, p2.y)
print(p3.name, p3.x, p3.y)

local x = 5

local p4 = p3 'p4' { x=\ x, y=\s 2*s.x }

print(p4.name, p4.x, p4.y)

local p5 = p4 'p5' { y={ 0, \s 3*s.x } }

print(p5.name, p5.x, p5.y[1], p5.y[2])

p3:set_function('f', \x,y (x+y, x-y))

print(p5.name, p5.f(2,3))

--p3:dump(nil, 4) -- 1st

print("+++ 2nd")
print(p3.name, p3.x)
p3.x = nil
print(p3.name, p3.x)

--p3:dump(nil, 4) -- 2nd

print("+++ 3rd")
rawset(p3, 'x', 3)   -- bypass [var]
print(p3.name, p3.x) -- no __index
p3.x = nil           -- true key delete
print(p3.name, p3.x) -- p3.__index -> p2.__index
p3.x = 3             -- true key delete
print(p3.name, p3.x) -- p3.__index

print(p5.toto)

--p3:dump(nil, 4) -- 3rd

-- TODO: test shallow copy...

--[=[
++ table: 0x0ea03638.name = nil
++ table: 0x0ea03638.x = 1
++ table: 0x0ea03638.y = 2
Point 1 2
++ p2.x = 2
++ p2.y = -1
p2  2 -1
++ p3.x = 3
++ p3.y = nil
++ p2.y = -1
p3  3 -1
++ p4.x = function: 0x0ea01898
++ p4.y = function: 0x0ea018d8
++ p4.x = function: 0x0ea01898
p4  5 10
++ p5.x = nil
++ p4.x = function: 0x0ea01898
++ p5.y = table: 0x0ea03e58
  ++ p5.table: 0x0ea02ab8 = nil
  ++ p4.table: 0x0ea02ab8 = nil
  ++ p3.table: 0x0ea02ab8 = nil
  ++ p2.table: 0x0ea02ab8 = nil
  ++ Point.table: 0x0ea02ab8 = nil
  ++ Object.table: 0x0ea02ab8 = nil

++ p5.y = table: 0x0ea03ec8
++ p5.x = nil
++ p4.x = function: 0x0ea01898

p5  5 0 15
++ p3.set_function = nil
++ p2.set_function = nil
++ Point.set_function = nil
++ Object.set_function = nil
++ p5.f = nil
++ p4.f = nil
++ p3.f = table: 0x0ea039e0
p5  5 -1
+++ 2nd
++ p3.x = 3
p3  3
++ p3.x = nil
++ p2.x = 2
p3  2
+++ 3rd
p3  3
++ p3.x = nil
++ p2.x = 2
p3  2
--]=]
