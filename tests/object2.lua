local Object = require 'object'

local Point = Object 'Point' {}

local p1 = Point      { x=1, y=2  }
local p2 = Point 'p2' { x=2, y=-1, z=0 }
local p3 = p2    'p3' { x=3  }

print(p1.name, p1.x, p1.y)
print(p2.name, p2.x, p2.y)
print(p3.name, p3.x, p3.y)

print(p3.name, p3:is_object(), p1.is_object(p3), p1.is_object({}))

-- example implementing notification

local function trace (fp, self, k, v)
  fp:write("object: '", self.name,
           "' is updated for key: '", tostring(k),
           "' with value: ")
  if type(v) == "string"
    then fp:write(": '", tostring(v), "'\n")
    else fp:write(":  ", tostring(v),  "\n") end
end

local function set_notification (self, file)
  local fp = file or io.stdout
  local nwidx = rawget(getmetatable(self) or {}, '__newindex')
  local mm = function (self, k, v)
    trace(fp, self, k, v) -- logging
    nwidx(    self, k, v) -- forward
  end
  self:set_metamethod('__newindex', mm, true) -- override!
end

set_notification(p2) -- new metamethod created, metatable is cloned

p2.x = 3 -- new behavior, notify about update 
p3.x = 4 -- created before set_metamethod (bad!), old behavior

local p4 = p2 'p4' { x=3 } -- new, inherit metatable
p4.x = 5 -- new behavior, notify about update

print(p2.name, p2.x, p2.y)
print(p3.name, p3.x, p3.y)
print(p4.name, p4.x, p4.y)
