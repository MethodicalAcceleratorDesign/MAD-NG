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
  The module object implements the necessary machinery to support prototype-
  based programming with some extensions:
  - On read, the lookup of values follows the inheritance down to 'Object'.
    If the retrieved value is a function, it is called with 'self' passed as
    argument (can be ignored) and the result is returned (i.e. property).
    To store functions with arguments and call semantic, use 'set_function'.
    To store metamethods for instances, use 'set_metamethod'.
  - On write, if the object has a defined __update metamethod, it will be
    called with 'self', key, and value as arguments, otherwise the value is
    stored.
  - __par points to the 'self' parent
  - __var points to the 'self' variables

RETURN VALUES
  The constructor of objects.

EXAMPLES
  Object = require 'object'
  Point = Object {}                     -- Point derives from Object
  p0 = Point { x=0, y=0 }               -- p0 is an instance of Point
  p1 = Point { x=1, y=1 }               -- p1 is an instance of Point
  p1.x, p1.y = 1, 2                     -- set p1.x and p1.y

SEE ALSO
  None
]]

-- implementation ------------------------------------------------------------o

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

-- metamethods
local meta = { -- from lj_obj.h
  '__add', '__call', '__concat', '__div', '__eq', '__gc', '__index', '__ipairs',
  '__le', '__len', '__lt', '__metatable', '__mod', '__mode', '__mul',
  '__newindex', '__pairs', '__pow', '__sub', '__tostring', '__unm', '__update',
}
for _,v in ipairs(meta) do meta[v]=v end -- build dictionary

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

local function is_notifiable (a)
  return rawget(getmetatable(a) or {}, '__update')
end

local function is_callable (a)
  return type(a) == 'function' or rawget(getmetatable(a) or {}, '__call')
end

-- objects are proxies controlling variables access and inheritance

function MT:__call (a) -- object ctor
  if is_string(a) then -- named obj
    local obj = {__id=a, [par]=self, [var]=var0, __index=self[var]}
    return setmetatable(obj, getmetatable(self))
  elseif is_table(a) then
    if self[var] == var0 then -- finalize named obj, set fast inheritance
      if rawget(self,'__id') ~= nil then a.__id, self.__id = self.__id, nil end
      self[var] = setmetatable(a, self);
      return self
    else -- unamed obj
      local obj = {[par]=self, [var]=a, __index=self[var]}
      setmetatable(a, obj) -- set fast inheritance
      return setmetatable(obj, getmetatable(self))
    end
  end
  error("invalid object argument")
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

function MT:__newindex (k, v) -- (+__update)
  local m = is_notifiable(self)
  if m then m(self, k, v) else self[var][k] = v end
end

function MT:__len   () return #self[var]        end
function MT:__pairs () return pairs(self[var])  end
function MT:__ipairs() return ipairs(self[var]) end

-- object methods

M.is_object = is_object

function M:isa (obj)
  assert(is_object(obj), "invalid argument")
  while self and self ~= obj do self = rawget(self,par) end
  return self == obj
end

function M:set_parent (obj)
  assert(self ~= M, "'Object' must be root")
  assert(is_object(obj), "invalid argument")
  self[par] = obj
  self.__index = obj[var]
  setmetatable(self, getmetatable(obj))
  return self
end

function M:set_function (name, func)
  assert(is_callable(func) or func == nil, "invalid argument")
  self[var][name] = is_function(func) and setmetatable({func}, MF) or func
  return self
end

function M:set_method (name, func)
  assert((is_callable(func) or func == nil) and meta[name] ~= name,
         "invalid argument")
  rawset(self, name, func)
  return self
end

function M:set_metamethod (name, func, override)
  assert((is_callable(func) or func == nil) and meta[name] == name
         and (MT[name] == nil or override), "invalid argument")
  local sm, pm = getmetatable(self), getmetatable(rawget(self,par))
  if sm == pm then -- create a new metatable if shared with parent
    sm={} ; for _,v in ipairs(meta) do sm[v] = pm[v] end
    setmetatable(self, sm)
  end
  rawset(sm, name, func)
  return self
end

--[=[ TODO
function M:make_class ()
  local sv = self[var] -- self[var] becomes self and discards old self
  sv[var] = setmetatable({},{__index=getmetatable(self)[var]});
  if rawget(sv,'__id') ~= nil then s[var].__id, sv.__id = sv.__id, nil end
  return setmetatable(sv, getmetatable(self)) -- set parent
end
--]=]

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

-- root object

M[var] = {
  __id ='Object',
  __par=\s rawget(s,par),
  __var=\s rawget(s,var),
  name =\s s.__id, -- alias
}

------------------------------------------------------------------------------o
return setmetatable(M,MT)

