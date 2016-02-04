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

local table  = require 'table.new'
local gmath  = require 'gmath'
local vector = require 'vector'

-- locals --------------------------------------------------------------------o

local max, floor = math.max, math.floor

-- implementation ------------------------------------------------------------o

function gmath.is_range (x)
  return getmetatable(x) == M
end

-- constructor

function range (start, stop, step)
  assert(start and step ~= 0, "invalid range argument")
  if not stop then start, stop = 0, start end
  if not step then step = start > stop and -1 or 1 end
  local size = max(0,floor((stop-start)/step+1.5))
  return setmetatable({_start=start, _stop=stop, _step=step, _size=size}, M)
end

-- methods

function M.range (r)
  return r._start, r._stop, r._step
end

function M.bounds (r)
  if r._step < 0 then
    return r._stop, r._start
  else
    return r._start, r._stop
  end    
end

function M.size (r)
  return r._size
end

function M.overlap (r, s)
  local rl, rh = r:bounds()
  local sl, sh = s:bounds()
  return not (rl < sl and rh < sl or rl > sh)
end

-- metamethods

M.__len = M.size

local function iter(r, i)
  if i < r._size then
    return i+1, r._start+i*r._step
  end
end

function M.__ipairs (r) -- iterator: for n in ipairs(r)
  return iter, r, 0
end

function M.__index (r, i)
  if type(i) == 'number' then
    return i >= 1 and i <= r._size and r._start+(i-1)*r._step or nil
  else
    return M[i]
  end
end

-- convertion

local function convert (r, ctor)
  local n = r._size
  local t = ctor(n,0)
  for i=1,n do
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

------------------------------------------------------------------------------o
return range
