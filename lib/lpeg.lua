local oss = ffi.os
local pcp = package.cpath
local plp = package.path
package.cpath = ";;./lib/lpeg/?-"..oss..".so;.\\lib\\lpeg\\?-"..oss..".dll;"
package.path = ";;./lib/lpeg/?.lua;.\\lib\\lpeg\\?.lua;"
local lpeg = require"lpeg"
local re   = require"re"
package.cpath = pcp
package.path  = plp

lpeg.setmaxstack(1024)

return lpeg, re
