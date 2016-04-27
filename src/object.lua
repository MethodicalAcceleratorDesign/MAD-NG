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
  The module object implements the necessary machinery to support prototype
  based programming with some extensions:
  - On read, value lookup follows the hierachy of objects, but functions are
    automatically called with 'self' passed as argument (can be ignored) before
    the value is returned. To avoid this behavior, use 'set_function' to store
    functions with arguments and call semantic.
  - On write, if value or new value have defined __copy metamethod, it will
    be called with both as arguments and the returned value will be stored.

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

-- special protected key to store object members
local var = setmetatable({}, {
  __index   =\ error "incomplete object initialization",
  __newindex=\ error "incomplete object initialization",
})

-- metatable of 'true' function proxy
local MF = {
  __call    =\t,... -> t[1](...),
  __index   =\ error "private object",
  __newindex=\ error "private const object",
}

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

local function is_copyable (a)
  return rawget(getmetatable(a) or {}, '__copy')
end

local function is_callable (a)
  return type(a) == 'function' or rawget(getmetatable(a) or {}, '__call')
end

-- objects are proxies controlling variables access and inheritance

local function obj2cls (self) -- convert object to a class
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

function M:__call (a) -- object ctor
  if is_string(a) then -- named obj
    return setmetatable({name=a, [var]=var}, obj2cls(self)) -- set parent
  elseif is_table(a) then
    if self[var] == var then -- finalize named obj, set inheritance
      if rawget(self,'name') ~= nil then a.name, self.name = self.name, nil end
      self[var] = setmetatable(a, {__index=getmetatable(self)[var]});
      return self
    else -- unamed obj
      setmetatable(a, {__index=self[var]}) -- set inheritance
      return setmetatable({[var]=a}, obj2cls(self)) -- set parent
    end
  end
  error("invalid object argument")
end

local function get (self, k) -- object lookup
  return self and (rawget(self,k) or get(getmetatable(self),k))
end

local function eval (self, v) -- variable lookup-eval
  return is_function(v) and v(self) or v
end

function M:__index (k) -- (+__eval+inheritance)
  return eval(self, self[var][k]) or get(getmetatable(self),k)
end

function M:__newindex (k, nv) -- (+__copy)
  local v = rawget(self[var],k)
  local c = v and (is_copyable(v) or is_copyable(nv))
  self[var][k] = c and c(v,nv) or nv
end

function M:__len    () return #self[var]        end
function M:__pairs  () return pairs(self[var])  end
function M:__ipairs () return ipairs(self[var]) end

-- object methods

function M:isa (obj)
  assert(is_object(obj), "invalid argument")
  while self and self ~= obj do
    self = getmetatable(self)
  end
  return self == obj
end

function M:set_function (name, func)
  assert(is_callable(func), "invalid argument")
  self[var][name] = is_function(func) and setmetatable({func}, MF) or func
  return self
end

function M:make_class ()
  local sv = self[var] -- self[var] becomes self and discards old self
  sv[var] = setmetatable({},{__index=getmetatable(self)[var]});
  if rawget(sv,'name') ~= nil then s[var].name, sv.name = sv.name, nil end
  return setmetatable(sv, getmetatable(self)) -- set parent
end

-- debug

function M:dump (file, level, indent, vars)
  local sv, pa, va = self[var], getmetatable(self), vars or {}
  local fp, lv, id = file or io.stdout, level or 1, indent or 1
  -- header
  for i=1,id-1 do fp:write("  ") end -- indent
  fp:write("objdump [", tostring(sv), "]\n")
  -- variables
  for k,v in pairs(sv) do
    for i=1,id do fp:write("  ") end -- indent
    fp:write(tostring(k))
        if is_string(v) then fp:write(": '", tostring(v)   , "'")
    elseif is_proxy (v) then fp:write(": [", tostring(v[1]), "]")
                        else fp:write(":  ", tostring(v)) end
    if va[k] and k ~= 'name' then fp:write(" (*)") else va[k] = true end
    fp:write("\n")
  end
  -- parent
  if lv > 1 and pa ~= M then pa:dump(fp, lv-1, id+1, va)
                        else fp:write("\n") end
end

-- root object

local object = setmetatable({
  [var] = { name  ='Object',
            parent=\s getmetatable(s),
            rawget=\s rawget(s,var),
          }}, M )

------------------------------------------------------------------------------o
return object
