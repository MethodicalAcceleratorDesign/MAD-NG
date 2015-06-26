package.path = '../?.lua;' .. package.path
local complex = require 'complex'

local a, b = (1+1i)/(math.sqrt(2)+1e-8), 1
local n = arg[1] and tonumber(arg[1]) or 1e8

for i=1,n do
	b = b * a
end

print('n=', n, 'a=', a, 'b=', b)
