-- time luajit -jv -Oloopunroll=50 -e "package.path = './?.lua;./lua/?.lua;./lib/?.lua;' .. package.path" lua/tests/matrix.lua

local gmath   = require 'gmath'
local complex = require 'complex'
local matrix  = require 'matrix'
local cmatrix = require 'cmatrix'
local jit    = require('jit')
local jitv   = require('jit.v')

jitv.start()

local sqrt, tostring = gmath.sqrt, gmath.tostring

local n = arg[1] and tonumber(arg[1]) or 1e7

local I, a, b, r
if true then
  I = 1.0/6
  a = matrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I} }
  b = matrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I} }
  r = matrix (6,6)
else
  I = complex(0,1.0/6)
  a = cmatrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I} }
  b = cmatrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I} }
  r = cmatrix (6,6)
end

-- check_bounds = true
-- print(a:get(100,100))

for i=1,n do
  b:mul(a,r)
  b, r = r, b
  -- b = b * a
end

print(tostring(b))
