local jit = require 'jit'
local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()
--jit.off()

local ffi = require 'ffi'
local complex = ffi.typeof 'complex'

local M={}

M.__new = function (ct, re, im)
  print("complex built")
  return setmetatable({re,im}, M)
end

M.__index = M

local function iter(r, i)
  if i < 1 then
    return i+1, r[i+1]
  end
end

function M.__ipairs (r) -- iterator: for n in ipairs(r)
  return iter, r, 0
end

-- local function complex(n)
--  return setmetatable({n=n, 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}, M)
-- end

-- needs to remove this line in the complex module...
ffi.metatype('complex', M)

local a = 0
for i=1,1e7 do
  local r = 1i
--  print(i,r, type(r))
  for _,v in ipairs(r) do
    a = a + v
  end
end
print(a)
