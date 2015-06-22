local pkg = "lpeg-0.12.2"
local oss = jit.os
local pcp = package.cpath
local plp = package.path
package.cpath = ";;./lib/"..pkg.."/?-"..oss..".so;.\\lib\\"..pkg.."\\?-"..oss..".dll;"
package.path = ";;./lib/"..pkg.."/?.lua;.\\lib\\"..pkg.."\\?.lua;"
local lpeg = require"lpeg"
local re   = require"re"
package.cpath = pcp
package.path  = plp

lpeg.setmaxstack(1024)

return lpeg, re
