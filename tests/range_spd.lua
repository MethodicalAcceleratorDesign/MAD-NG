local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()

local range = require 'range'

local a = 0
local r = range(1e9)
print(r:range())

for _,v in ipairs(r) do
  a = a + v
end

print(a)
