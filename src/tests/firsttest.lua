luaunit = require('luaunit')

--[[
function add(v1,v2)
    -- add positive numbers
    -- return 0 if any of the numbers are 0
    -- error if any of the two numbers are negative
    if v1 < 0 or v2 < 0 then
        error('Can only add positive or null numbers, received '..v1..' and '..v2)
    end
    if v1 == 0 or v2 == 0 then
        return 0
    end
    return v1+v2
end
--]]

local gmath = require "gmath"

local add = gmath.add

function testAddPositive()
    luaunit.assertEquals(add(1,1),2)
end

function testAddZero()
    luaunit.assertEquals(add(1,0),1)
    luaunit.assertEquals(add(0,5),5)
    luaunit.assertEquals(add(0,0),0)
end

local complex = require "complex"
local sqrt = gmath.sqrt

function testComplexSqrt()
    luaunit.assertEquals(sqrt(complex(-1,0)),1i)
end

os.exit( luaunit.LuaUnit.run() )