--[[
Wrapper for lpeg.so.
Requires and returns lpeg.so, keeping the cpath as it is.
]]

local pcp = package.cpath
package.cpath = ";;./lib/lpeg/lpeg-0.12/?.so;.\\lib\\?\\lpeg-0.12\\?.dll;"
local lpeg = require"lpeg"
lpeg.setmaxstack(1024)
package.cpath = pcp
return lpeg
