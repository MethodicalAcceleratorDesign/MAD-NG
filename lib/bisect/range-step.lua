-- bad commit ccae333844c7aad0934f13f7698894c883a6b561

local ffi = require 'ffi'
ffi.cdef[[
  typedef struct range { double start, stop, step; } range;
]]
local range = ffi.typeof 'range'
local a = 1
local b = 0
for i=1,100 do
  r = range(0,-1,-a)
  b = r.start + r.step*1
--  print( i, a, r.start, r.step, r.start + r.step*1 )
end

--print(b)
os.exit(b == -1)