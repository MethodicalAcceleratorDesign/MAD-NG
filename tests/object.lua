local Object = require 'object'

local Point = Object 'Point' {}

local p1 = Point      { x=1, y=2  }
local p2 = Point 'p2' { x=2, y=-1, z=0 }
local p3 = p2    'p3' { x=3  }

print(p1.name, p1.x, p1.y)
print(p2.name, p2.x, p2.y)
print(p3.name, p3.x, p3.y)

local x = 5

local p4 = p3 'p4' { x=\s s.parent.x, y=\s 2*s.x }

print(p4.name, p4.x, p4.y, p4.z)
print(p4.name, p4.rawget.x, p4.rawget.y, p4.rawget.z)

local p5 = p4 'p5' { y=\ { 0, 3*x } }

print(p5.name, p5.x, p5.y[1], p5.y[2])

local p6 = p5 'p6' { x=5 }

print(p6.name, p6.x, p6.y[1], p6.y[2])

p3:set_function('f', \x,y (x+y, x-y))

print(p5.name, p5.f(2,3))

p3:dump(nil, 4) -- 1st

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

print("list of p5 variables")
for k,v in pairs(p5) do
  print(k, ":", v, p5[k])
end 

--p3:dump(nil, 4) -- 3rd
