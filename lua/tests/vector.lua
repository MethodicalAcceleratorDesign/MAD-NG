local generic = require 'generic'
local complex = require 'complex'
local vector  = require 'vector'

local sqrt = generic.sqrt

local n = arg[1] and tonumber(arg[1]) or 1e7

local I = 1 -- complex(0,1)
local a = vector {1.0001, 1.0002, 1.0003, 1.0004}
local b = vector {I,I,I,I}

for i=1,n do
	b = a * b + b - a
end

for i=1,#b do io.write(generic.tostring(b[i]), ' ') end
io.write('\n')

local x = vector {1,0,0}
local y = vector {0,1,0}

print( x:angle(y), y:angle(x) )
