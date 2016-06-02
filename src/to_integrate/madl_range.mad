--[=[
 o----------------------------------------------------------------------------o
 |
 | Range module
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
  
  Purpose:
  - provides full set of functions and operations on ranges

  Information:
  - ranges are not tables nor vector. totable and tovector convert a range in
    a table or a vector.

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  range

SYNOPSIS
  local range = require 'range'
  TODO

DESCRIPTION
  Ranges describe start, stop (included) and step.

REMARK:
  - Indexing a range outside its bounds return nil.
  - Ranges can be used as (stateless) iterators in for loops 

RETURN VALUES
  The constructor of ranges

SEE ALSO
  vector, cvector, matrix, cmatrix
]]
 
-- modules -------------------------------------------------------------------o

local ffi    = require 'ffi'
local gmath  = require 'gmath'
local vector = require 'vector'
local table  = require 'table.new'

-- ffi -----------------------------------------------------------------------o

ffi.cdef[[
typedef struct { double _start, _stop, _step, _size; } range_t;
]]

-- locals --------------------------------------------------------------------o

local istype = ffi.istype
local max, floor, isnum = math.max, math.floor, gmath.is_number

-- FFI type constructors
local range_ctor = ffi.typeof 'range_t'

-- implementation ------------------------------------------------------------o

function gmath.is_range (x)
  return type(x) == 'cdata' and istype('range_t', x)
end

-- constructor

local function size (start, stop, step)
--  if step < -1 or step > 1 then
--    return max(0, floor( (stop-start)/step +step +0.5) )
--  end
  if step > 1 then
    return max(0, floor( (stop-start+1)/step +0.5) )
  elseif step < -1 then
    return max(0, floor( (stop-start-1)/step +0.5) )
  else
    return max(0, floor( (stop-start  )/step +1.5) )
  end
end
  
function range (start, stop, step)
  assert(start and step ~= 0, "invalid range argument")
  local r = range_ctor()
  if not stop then start, stop = 0, start end
  if not step then step = start > stop and -1 or 1 end
  local size = size(start, stop, step)
  r._start, r._stop, r._step, r._size = start, stop, step, size
  return r
end

-- mutator

function M.scale (r, a)
  r._stop, r._step = r._stop*a, r._step*a
  r._size = size(r._start, r._stop, r._step)
  return r
end

function M.translate (r, a)
  r._start, r._stop = r._start+a, r._stop+a
  r._size = size(r._start, r._stop, r._step)
  return r
end

function M.reverse (r)
  r._start, r._stop, r._step = r._stop, r._start, -r._step
  return r
end

-- methods

function M.range (r)
  return r._start, r._stop, r._step
end

function M.size (r)
  return r._size
end

function M.first (r)
  return r._size > 0 and r._start or nil
end

function M.last (r)
  return r:value(r._size)
end

--[[function M.index (r, x)
  local i = floor((x-r._start)/r._step+1.5)
  return i >= 1 and i <= r._size and i or nil
end]]

function M.index (r, x)
  local i = size(r._start, x, r._step)
  return i >= 1 and i <= r._size and i or nil
end

function M.value (r, i)
  return i >= 1 and i <= r._size and r._start+(i-1)*r._step or nil
end

function M.element (r, x)
  return x == r:value(r:index(x))
end

function M.minmax (r)
  if r._step < 0 then
    return r:last(), r:first()
  else
    return r:first(), r:last()
  end    
end

function M.bounds (r)
  if r._step < 0 then
    return r._stop, r._start
  else
    return r._start, r._stop
  end    
end

function M.overlap (r, s)
  if r._size == 0 or s._size == 0 then
    return false
  end

  local rl, rh = r:bounds()
  local sl, sh = s:bounds()
  return not (rl < sl and rh < sl or rl > sh)
end

-- convertion

local function convert (r, ctor)
  local t = ctor(r._size,0)
  for i=1,r._size do
    t[i] = r._start+(i-1)*r._step
  end
  return t
end

function M.totable (r)
  return convert(r, table)
end

function M.tovector (r)
  return convert(r, vector)
end

function M.tostring (r)
  if r._step == 1 then
    return string.format("%g:%g", r._start, r._stop)
  else
    return string.format("%g:%g:%g", r._start, r._stop, r._step)
  end
end

function M.equal (r1, r2)
  if r1._size == 0 and r2._size == 0 then
    return true
  else
    return r1._start == r2._start and r1._step == r2._step and r1._size == r2._size
  end          
end

-- metamethods

M.__len      = M.size
M.__tostring = M.tostring
M.__eq       = M.equal

local function iter(r, i)
  if i < r._size then
    return i+1, r._start+i*r._step
  end
end

function M.__ipairs (r) -- iterator: for n in ipairs(r)
  return iter, r, 0
end

function M.__index (r, i)
  return isnum(i) and r:value(i) or M[i]
end




ffi.metatype('range_t', M)

------------------------------------------------------------------------------o
return range
