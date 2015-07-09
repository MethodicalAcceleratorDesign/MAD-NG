local generic = require 'generic'
local complex = require 'complex'
local matrix  = require 'matrix'

local sqrt = generic.sqrt

local n = arg[1] and tonumber(arg[1]) or 1e7

local I = 1 -- complex(0,1)
local a = matrix { {1.0001, 1.0002}, {1.0003, 1.0004} }
local b = matrix { {I,I}, {I,I} }

for i=1,n do
	b = a * b + b - a
end

print(b:tostring())
