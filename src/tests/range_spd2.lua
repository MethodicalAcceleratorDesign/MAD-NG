local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()

local range = require 'range'

local a = 0
for i=1,1e8 do
  local r = range(10)
  for _,v in ipairs(r) do
    a = a + v
  end
end
print(a)
