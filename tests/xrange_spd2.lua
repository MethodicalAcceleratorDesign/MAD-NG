local jit = require 'jit'
local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()
--jit.off()

local range = require 'xrange'

local a = 0
for i=1,100 do
  local r = 0..2
  a = a + r:len()
  print(i, r)
  print(r.start, r.stop, r.step)
end
print(a)
