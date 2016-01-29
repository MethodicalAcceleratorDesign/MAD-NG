--[=[
 o----------------------------------------------------------------------------o
 |
 | Object module (object model)
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
  - TODO

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- module --------------------------------------------------------------------o

M.__help.self = [[
NAME
  object -- creates objects

SYNOPSIS
  object = require"object"
  obj1 = object {}               -- create a new empty object
  obj2 = object { ... }          -- create a new object with values

DESCRIPTION
  The module object creates new objects from lists that become instances
  of their parent, with callable semantic (i.e. constructor). A 'list' is a
  table without metatable.
  
  The returned object has its parent (its constructor) set as metatable and
  inherits all properties of it automatically, hence implementing a prototype
  language.

  Hence an object is a 'table' that can be used as a constructor (a function)
  to create new instances of itself, an object (a table) or a class/parent
  (a metatable).

  The module provides some utilities like :spr, :isa, and .is_object for type
  identification, :set and :cpy for data manipulation, and .is_list to check
  for lists, that is 'virgin' tables.

RETURN VALUES
  The input list properly setup to be an object.

ERRORS
  If the object does not receive a list, an invalid argument error is raised.

EXAMPLES
  Object = require 'object'
  Point = Object {}                     -- Point derives from Object
  p0 = Point { x=0, y=0 }               -- p0 is an instance of Point
  p1 = Point { x=1, y=1 }               -- p1 is an instance of Point
  p2 = p1:cpy()                         -- p2 is a copy of p1
  p1:set { x=-1, y=-2 }                 -- set p1.x and p1.y (slow)
  p1.x, p1.y = 1, 2                     -- set p1.x and p1.y (faster)

  is_list = require"utils".is_list
  is_list { x=0, y=0 }                  -- return true
  is_list (p0)                          -- return false
  p0:isa(Point)                         -- return true

SEE ALSO
  None
]]

-- locals ----------------------------------------------------------------------

local getmetatable, setmetatable = getmetatable, setmetatable
local type, pairs = type, pairs
local is_list = require "utils" .is_list

local MT = {}; setmetatable(M, MT) -- make this module the root of all objects

-- members ---------------------------------------------------------------------

M.name = 'object'
M.is_object = true

-- methods ---------------------------------------------------------------------

-- return the next parent
function M:spr()
  return getmetatable(self)
end

-- return the parent id or nil
function M:isa(id)
  local a = getmetatable(self);
  while a ~= nil and a ~= id do a = getmetatable(a) end
  return a ~= nil
end

-- set values taken from iterator
function M:set(a)
  for k,v in pairs(a) do self[k] = v end
  return self
end

-- make a copy
function M:cpy()
  return setmetatable({}, getmetatable(self)):set(self)
end

-- metamethods -----------------------------------------------------------------

-- constructor
function MT:__call(a)
  if is_list(a) then
    if not rawget(self, '__call') then
      self.__index = self         -- inheritance
      self.__call  = MT.__call    -- constructor
    end
    return setmetatable(a, self)
  end
  error ("invalid constructor argument, list expected")
end

-- end -------------------------------------------------------------------------
return M
