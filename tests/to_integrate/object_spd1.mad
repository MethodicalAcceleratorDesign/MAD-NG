local Object = require 'object'

local Point = Object 'Point' {}

local p1 = Point 'p1' { x=3, y=2, z=1  }
local p2 = p1    'p2' { x=2, y=1 }
local p3 = p2    'p3' { x=1  }
local p4 = p3    'p4' { }

print(p1.name, p1.x, p1.y, p1.z)
print(p2.name, p2.x, p2.y, p2.z)
print(p3.name, p3.x, p3.y, p3.z)
print(p4.name, p4.x, p4.y, p4.z)

local s = 0
for i=1,5e8 do
  s = s + (p1.z + p2.z + p3.z + p4.z)
end

print('s=',s)
