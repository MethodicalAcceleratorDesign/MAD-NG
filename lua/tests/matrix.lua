-- time luajit -jv -Oloopunroll=50 -e "package.path = './?.lua;./lua/?.lua;./lib/?.lua;' .. package.path" lua/tests/matrix.lua

local complex = require 'complex'
local matrix  = require 'matrix'
local cmatrix = require 'cmatrix'
local gmath   = require 'gmath'

local sqrt, tostring = gmath.sqrt, gmath.tostring

local n = arg[1] and tonumber(arg[1]) or 1e7
--1.0/6 -- 
local I = complex(0,1.0/6)
local a = cmatrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I},}
local b = cmatrix { {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I}, {I,I,I,I,I,I},}
local r = cmatrix (6,6)

for i=1,n do
  b:mul(a,r)
  b, r = r, b
  -- b = b * a
end

print(tostring(b))
