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
  - Provide an object model to support prototype-based programming with value
    semantic for functions stored in variables and further extensions. 

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
  The 'object' module implements the necessary machinery to support prototype-
  based programming with value semantic for functions and further extensions:
  - On read, the lookup of values follows the inheritance down to 'Object' with
    precedence of variables over methods.
    + If the retrieved value is a function, it is called with 'self' passed as
      argument (can be ignored) and it returns the result.
    + To store functions with arguments in variables, use 'set_function'.
    + To store methods or metamethods for instances (both are inherited), use
      'set_method' or 'set_metamethod' respectively.
  - On write, the value is simply stored (no lookup).
    + To override this behavior, just (re)defined the __newindex metamethod
      using set_metamethod with 'override' as true (use with care!).
  - On build, the new instance is connected to its parent (inheritance).
    + If the new instance has a defined __init metamethod (inherited), it will
      be called on the new instance and non-nil result is returned.
  - Root 'Object' defines the following variables:
    + 'name'  points to 'self' name unless overridden
    + '__par' points to 'self' parent unless overridden
    + '__var' points to 'self' variables unless overridden
    + '__id'  holds 'self' name (may be inherited)

RETURN VALUES
  The constructor of objects.

EXAMPLES
  Object = require 'object'
  Point = Object {}              -- point is an instance of Object
  p0 = Point { x=0, y=0 }        -- p0 is an instance of Point
  p1 = Point { x=1, y=1 }        -- p1 is an instance of Point
  p2 = p1 { x=3 }                -- p2 is an instance of p1 and inherits p1.y
  p1.x, p1.z = 2, 3              -- set p1.x, p1.z
  print(p1.x, p1.y, p1.z)        -- print 2 1 3 
  print(p2.x, p2.y, p2.z)        -- print 3 1 3 

SEE ALSO
  None
]]

-- documentation -------------------------------------------------------------o

--[=[
  Schematic object-model representation:
  --------------------------------------

  o0 = require 'object'
  o1 = o0 'obj1' {*o1-var*}
  o2 = o1 'obj2' {*o2-var*}
  o3 = o1 'obj2' {*o3-var*}                    +-------------+
                             +---------------+>| *meta-tbl*  |<------------+
+---------+                  |   +---------+ | | metamethods | +---------+ |
|  *o2*   |                  |   |  *o1*   | | +-------------+ |  *o0*   | |
|  [meta] |------------------+   |  [meta] |-+                 |  [meta] |-+
|   [par] |------------------|-->|   [par] |------------------>|   [par] |-->.
| __index |------------------|-+ | __index |-----------------+ | __index |-->.
|         |  +-----------+   | | |         |  +-----------+  | |         |
|   [var] |->|  *o2-var* |   | | |   [var] |->| *o1-var*  |  | |   [var] |-+
| methods |  |    [meta] |-+ | | | methods |  |    [meta] |-+| | methods | |
+---------+  | variables | | | | +---------+  | variables | || +---------+ |
     ^       +-----------+ | | |      ^       +-----------+ ||        +----+
     +---------------------+ | |      |             ^       ||        v
+---------+                  | |      +-------------|-------+|  +-----------+
|  *o3*   |                  | |      |             |        |  | *o0-var*  |
|  [meta] |------------------+ |      |             |        +->|    [meta] |->.
|   [par] |--------------------|------+             |           | variables |
| __index |--------------------+--------------------+           +-----------+
|         |  +-----------+
|   [var] |->| *o3-var*  |
| methods |  |    [meta] |-+
+---------+  | variables | |
     ^       +-----------+ |
     +---------------------+

  Catching creation:
  ------------------

  Example how to count the number of objects created

  local count = 0
  local function set_counter (self)
    local mm = function (self)
      count = count + 1
      return self
    end
    self:set_metamethod('__init', mm)
  end
  set_counter(o0) -- before o1,o2,o3 creation

  Catching writes:
  ----------------

  Example how to set a notification-on-write with logging

  local function set_notification (self)
    local nwidx = rawget(getmetatable(self),'__newindex')
    local mm = function (self, k, v)
      trace(self, k, v) -- logging
      nwidx(self, k, v) -- forward
    end
    self:set_metamethod('__newindex', mm, true) -- override!
  end
  set_notification(o1) -- before o2,o3 creation
]=]

-- implementation ------------------------------------------------------------o

-- metamethods
local meta = { -- from lj_obj.h
  '__add', '__call', '__concat', '__div', '__eq', '__gc', '__index', '__init',
  '__ipairs', '__le', '__len', '__lt', '__metatable', '__mod', '__mode',
  '__mul', '__new', '__newindex', '__pairs', '__pow', '__sub', '__tostring',
  '__unm',
}
for _,v in ipairs(meta) do meta[v]=v end -- build dictionary

-- special protected key to store parent and object members
local par, var, var0 = {}, {}, setmetatable({}, {
  __index   =\ error "incomplete object initialization",
  __newindex=\ error "incomplete object initialization",
})

-- metatable of 'true' function proxy
local MF = {
  __call    =\t,... -> t[1](...),
  __index   =\ error "private object",
  __newindex=\ error "private const object",
}

-- metatable of 'Object'
local MT = {}

-- helpers

local function is_string (a)
  return type(a) == 'string'
end

local function is_function (a)
  return type(a) == 'function'
end

local function is_proxy (a)
  return type(a) == 'table' and getmetatable(a) == MF
end

local function is_table (a)
  return type(a) == 'table' and getmetatable(a) == nil
end

local function is_object (a)
  return type(a) == 'table' and getmetatable(a) ~= nil and rawget(a,var) ~= nil
end

local function is_callable (a)
  return type(a) == 'function' or rawget(getmetatable(a) or {}, '__call')
end

-- objects are proxies controlling variables access and inheritance

local function init(a)
  local m = rawget(getmetatable(a), '__init')
  return m and m(a) or a
end

function MT:__call (a) -- object ctor
  if is_string(a) then -- named obj
    local obj = {__id=a, [par]=self, [var]=var0, __index=self[var]}
    return setmetatable(obj, getmetatable(self)) -- incomplete obj
  elseif is_table(a) then
    if self[var] == var0 then -- finalize named obj
      if rawget(self,'__id') ~= nil then a.__id, self.__id = self.__id, nil end
      self[var] = setmetatable(a, self); -- set fast inheritance
      return init(self)
    else -- unnamed obj
      local obj = {[par]=self, [var]=a, __index=self[var]}
      setmetatable(a, obj) -- set fast inheritance
      return init(setmetatable(obj, getmetatable(self))) -- complete obj
    end
  end
  error("invalid object constructor argument, string or table expected")
end

local function eval (self, v) -- variable eval
  return is_function(v) and v(self) or v
end

local function get (self, k) -- object lookup
--  if self == nil then return nil else
--  io.write("get  : ", k, " '", tostring(self[var].__id or 'nil'), "'\n") end
  return self and (rawget(self,k) or get(rawget(self,par),k))
end

function MT:__index (k) -- (+eval+inheritance)
--  io.write("index: ", k, " '", tostring(self[var].__id or 'nil'), "'\n")
  return eval(self, self[var][k]) or get(rawget(self,par),k)
end

function MT:__newindex (k, v)
  self[var][k] = v
end

function MT:__len   () return #self[var]        end
function MT:__pairs () return pairs(self[var])  end
function MT:__ipairs() return ipairs(self[var]) end

-- object methods

M.is_object = is_object

function M:isa (obj)
  assert(is_object(obj), "invalid 'obj' argument, valid object expected")
  while self and self ~= obj do self = rawget(self,par) end
  return self == obj
end

function M:set_function (name, func)
  assert(is_object(self)               , "invalid 'self' argument, valid object expected")
  assert(is_callable(func) or func==nil, "invalid 'func' argument, not callable")
  self[var][name] = is_function(func) and setmetatable({func}, MF) or func
  return self
end

function M:set_method (name, func)
  assert(is_object(self)               , "invalid 'self' argument, valid object expected")
  assert(meta[name] ~= name            , "invalid 'name' argument, metamethod detected")
  assert(is_callable(func) or func==nil, "invalid 'func' argument, not callable")
  rawset(self, name, func)
  return self
end

function M:set_metamethod (name, func, override)
  assert(is_object(self)               , "invalid 'self' argument, valid object expected")
  assert(meta[name] == name            , "invalid 'name' argument, not a metamethod")
  assert(is_callable(func) or func==nil, "invalid 'func' argument, not callable")
  local sm, pm = getmetatable(self), getmetatable(rawget(self,par))
  assert(rawget(sm, name) == nil or override, "cannot override inherited behavior")
  if sm == pm then -- create a new metatable if shared with parent
    sm={} ; for _,v in ipairs(meta) do sm[v] = pm[v] end
    setmetatable(self, sm)
  end
  rawset(sm, name, func)
  return self
end

function M:set_parent (obj, name)
  assert(self ~= M                        , "'Object' must stay the root of objects")
  assert(is_object(self) or is_table(self), "invalid 'self' argument, table or object expected")
  assert(is_object(obj)                   , "invalid 'obj' argument, valid object expected")
  rawset(self, par, obj)
  rawset(self, '__index', obj[var])
  if is_table(self) and rawget(self,var) == nil then
    rawset(self, var, setmetatable({__id=name},self))
  end 
  setmetatable(self, getmetatable(obj))
  return self
end

-- debug

function M:dump (file, level, indent, vars)
  local sv = self[var] 
  file, level, indent, vars =
    file or io.stdout, level or 1e6, indent or 1, vars or {}
  -- header
  for i=1,indent-1 do file:write("  ") end -- indent
  file:write("+ [", tostring(sv), "]\n")
  -- variables
  for k,v in pairs(sv) do
    for i=1,indent do file:write("  ") end -- indent
    file:write(tostring(k))
        if is_string(v) then file:write(": '", tostring(v)   , "'")
    elseif is_proxy (v) then file:write(": [", tostring(v[1]), "]")
                        else file:write(":  ", tostring(v)) end
    if vars[k] and string.sub(k,1,2) ~= '__' then
      file:write(" (")
      for i=1,vars[k] do file:write('*') end
      file:write(")");
    end
    vars[k] = (vars[k] or 0)+1
    file:write("\n")
  end
  -- parent
  if level > 1 and self ~= M then
    self[par]:dump(file, level-1, indent+1, vars)
  end
end

-- root Object = module

M[var] = {
  __id ='Object',
  __par=\s rawget(s,par),
  __var=\s rawget(s,var),
  name =\s s.__id, -- alias
}

------------------------------------------------------------------------------o
return setmetatable(M,MT)

