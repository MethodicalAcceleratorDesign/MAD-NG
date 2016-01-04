local complex = require 'complex'
local matrix  = require 'matrix'
local gmath   = require 'gmath'

local sqrt, tostring = gmath.sqrt, gmath.tostring

local n = arg[1] and tonumber(arg[1]) or 1e7

local function vector(a) return matrix {a} end

local I = 1 -- complex(0,1)
local a = vector {1.0001, 1.0002, 1.0003, 1.0004}
local b = vector {I,I,I,I}

for i=1,n do
	b = a * b + b - a
end

io.write(tostring(b), '\n')

local x = vector {1,0,0}
local y = vector {0,1,0}

print( x:angle(y), y:angle(x) )
