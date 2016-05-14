local range = require 'xrange'

local a = 0

local t = tonumber(arg[1]) or 1
io.write('test: ', t, '\n')

if t == 1 then
-- 2.4 sec
for i,v in ipairs(0..1e9) do
  a = a + i - v
end

elseif t == 2 then
-- 2.4 sec
for i,v in ipairs(0..9e8..0.9) do
  a = a + i - v
end

elseif t == 3 then
-- 2.8 sec
local r = 0..1e9
for i=1,1e9 do
  a = a + i - r[i]
end

elseif t == 4 then
-- 3.9 sec
for k=1,1e7 do
 for i,v in ipairs(0..100) do
   a = a + i - v
 end
end

end

io.write('sum : ', a, '\n')
