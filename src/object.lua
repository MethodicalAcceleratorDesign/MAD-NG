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
  - Provide an object model to support prototype-based programming

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  object -- creates objects

SYNOPSIS
  object = require 'object'
  obj1 = object {}               -- create a new empty object
  obj2 = object { ... }          -- create a new object with values
  obj3 = object 'name' { ... }   -- create a new object with name and values

DESCRIPTION
  The module object implements the necessary machinery to support
  prototype-based programming and deferred expressions.
  
RETURN VALUES
  The constructor of objects.

EXAMPLES
  Object = require 'object'
  Point = Object {}                     -- Point derives from Object
  p0 = Point { x=0, y=0 }               -- p0 is an instance of Point
  p1 = Point { x=1, y=1 }               -- p1 is an instance of Point
  p1:set { x=-1, y=-2 }                 -- set p1.x and p1.y (slow)
  p1.x, p1.y = 1, 2                     -- set p1.x and p1.y (faster)

SEE ALSO
  None
]]

-- implementation ------------------------------------------------------------o

local MT  = {} -- metatable of object (root)
local MTT = {} -- metatable of tables (proxy)
 
local var = {} -- special key to store object members
local obj = {} -- special key to store object self reference

-- protect var

setmetatable(var, {
  __index    =\ error "incomplete object initialization",
  __newindex =\ error "incomplete object initialization"
})

-- helpers

local function is_string (a)
  return type(a) == 'string'
end

local function is_function (a)
  return type(a) == 'function'
end

local function is_table (a)
  return type(a) == 'table' and getmetatable(a) == nil
end

local function init_parent (self) -- init parent as a class
  if rawget(self, '__call') == nil then
    local mt = getmetatable(self) -- copy metamethods
    rawset(self, '__call'    , rawget(mt, '__call'    ))
    rawset(self, '__index'   , rawget(mt, '__index'   ))
    rawset(self, '__newindex', rawget(mt, '__newindex'))
    rawset(self, '__len'     , rawget(mt, '__len'     ))
    rawset(self, '__pairs'   , rawget(mt, '__pairs'   ))
    rawset(self, '__ipairs'  , rawget(mt, '__ipairs'  ))
  end
  return self
end

local function init_class (self) -- init object as a class
  local sv = self[var]
  sv[var] = {}                                -- set var
  if rawget(self, name) ~= nil then
    sv.name = self.name                       -- set name
  end
  return setmetatable(sv, getmetatable(self)) -- set parent
end

local function init_table (self, k, v) -- put table in a proxy
  local sv = self[var]
  sv[k] = setmetatable({[var]=v, [obj]=self[obj] or self}, MTT)
  return sv[k]
end

local function eval_tbl(self, k, v) -- read-eval table members
  return is_function(v) and v(self[obj]) or
         is_table   (v) and init_table(self, k, v) or v
end

local function eval_obj(self, k, v) -- read-eval object members
  return is_function(v) and v(self) or
         is_table   (v) and init_table(self, k, v) or v
end

-- proxy for controlling tables stored in object variables

function MTT:__index (k)
  local v = self[var][k]
  return v and eval_tbl(self, k, v)
end

function MTT:__newindex (k, v)   self[var][k] = v  end
function MTT:__len    () return #self[var]        end
function MTT:__pairs  () return pairs(self[var])  end
function MTT:__ipairs () return ipairs(self[var]) end

-- proxy for controlling object variables

function MT:__call (a) -- object ctor
  if is_table(a) then
    if self[var] == var then self[var] = a; return self end
    return setmetatable( {        [var]=a  }, init_parent(self) )
  elseif is_string(a) then
    return setmetatable( {name=a, [var]=var}, init_parent(self) )
  end
  error("invalid object argument")
end

function MT:__index (k) -- (+inheritance)
  local v = self[var][k]
  return v and eval_obj(self, k, v) or getmetatable(self)[k]
end

function MT:__newindex (k, v)   self[var][k] = v  end
function MT:__len    () return #self[var]         end
function MT:__pairs  () return pairs(self[var])   end
function MT:__ipairs () return ipairs(self[var])  end
function MT:parent   () return getmetatable(self) end

function MT:get (k) -- idem __index (-eval)
  return self[var][k] or getmetatable(self)[k]
end

function MT:set (a, v) -- idem __newindex (+shallow copy)
  local sv = self[var]
  if v ~= nil then
    sv[a] = v
  else -- shallow copy (slow)
    for k,v in pairs(a) do sv[k] = v end
  end
  return self
end

function MT:set_method(k, f)
  if is_function(f) then
    rawset(self,k,f)
    return self
  end
  error("invalid set_method argument")
end

function MT:make_class ()
  if self[var] ~= nil then
    return self:init_class()
  end
  error("invalid or incomplete object")
end

function MT:dump (file, level, indent)
  local fp = file or io.stdout
  local lv = level or 1
  local id = indent or 1
  local sv = self[var]
  local pa = self:parent()
  if id == 1 then
    fp:write("objdump '", self.name, "' [", tostring(self), "]\n")
  end
  for k,v in pairs(sv) do
    for i=1,id do fp:write("  ") end
    fp:write(tostring(k), ": ", tostring(v), "\n")
  end
  if lv > 1 and pa ~= MT then
    for i=1,id do fp:write("  ") end
    fp:write("parent '", pa.name, "' [", tostring(pa), "]\n")
    pa:dump(fp, lv-1, id+1)
  end
  fp:write("\n")
end

------------------------------------------------------------------------------o
return setmetatable( {name='Object', [var]={}}, MT ) -- root
