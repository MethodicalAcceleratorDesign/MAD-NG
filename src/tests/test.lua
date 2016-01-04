local ffi = require 'ffi'

ffi.cdef[[
	void myprint(void);
	void mycount(void);
]]

local clib = ffi.C

print(math.cos(math.atan(1)*2))
