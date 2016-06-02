local Object = require 'object'

-- example implementing notification

local count = 0
local function set_counter (self)
  local mm = function (self)
    count = count + 1
    return self
  end
  self:set_metamethod('__init', mm)
end
set_counter(Object)

local Point = Object 'Point' {}

local p1 = Point      { x=1, y=2  }
local p2 = Point 'p2' { x=2, y=-1, z=0 }
local p3 = p2    'p3' { x=3  }

print(p1.name, p1.x, p1.y)
print(p2.name, p2.x, p2.y)
print(p3.name, p3.x, p3.y)
print(p3.name, p3:is_object(), p1.is_object(p3), p1.is_object({}))

print(count)