package.path = '../?.lua;' .. package.path
local complex = require 'complex'

local a, b = (1+1i)/(math.sqrt(2)+1e-8), 1
local n = arg[1] and tonumber(arg[1]) or 1e8

for i=1,n do
	b = b * a
end

print('n=', n, 'a=', a, 'b=', b)
print('a^2=', a^2)
print('a^n=', a^n)
print('a^0.5=', a^0.5)
print('sqrt(a)=', a:sqrt())

--[[
mathematica:
0.493068692[627911512064026273807993328818603327910066849118793 + 0i
online accurate calc:
0.493068692[627911512064026273807993328818603 + 0i
hp42s:
0.493068692[628    +  0i

brute force (1e8 loops):
0.493068687[91051  +  0i
cpow (direct ^1e8):
0.493068687[15417  -  4.8176420190897e-09i
ipow (log2(1e8) ~ 27 loops):
0.493068687[57683  +  0i
]]

