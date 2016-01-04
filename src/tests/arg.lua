ffi = require "ffi"

print(...)

for k,v in pairs(arg) do
	io.write('arg[', k, ']=', v, '\n')
end

ffi.cdef[[
void mad_fatal(const char*);
void mad_error(const char*);
void mad_lua_setloc(int level);
]]

local function myerror()
--	ffi.C.mad_lua_setloc(2)
--	ffi.C.mad_error("This is an error")

	ffi.C.mad_lua_setloc(2)
	ffi.C.mad_fatal("This is a fatal error")
end

myerror()
-- print('returned from mad_error')

error("This is a fatal error")

ffi.C.mad_lua_setloc(1); ffi.C.mad_error("This is a fatal error")