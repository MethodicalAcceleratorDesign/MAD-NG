-- local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
-- jitv.start()

local range = require 'xrange'

local a = 0

-- 2.4 sec
for i,v in ipairs(0..1e9) do
  a = a + i - v
end

-- 2.4 sec
-- for i,v in ipairs(0..9e8..0.9) do
--   a = a + i - v
-- end

-- 2.8 sec
-- local r = 0..1e9
-- for i=1,1e9 do
--   a = a + i - r[i]
-- end

-- 3.9 sec
-- for k=1,1e7 do
--  for i,v in ipairs(0..100) do
--    a = a + i - v
-- -- print(i, v)
--  end
-- end

print(a)
