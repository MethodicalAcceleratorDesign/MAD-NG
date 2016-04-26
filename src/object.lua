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
  prototype-based programming, and __eval and __copy metamethods.
  
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

local var = {} -- special key to store object members
local obj = {} -- special key to store object self reference
local MT  = {} -- metatable of tables (proxy)
 
-- protect var

setmetatable(var, {
  __index    =\ error "incomplete object initialization",
  __newindex =\ error "incomplete object initialization",
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

local function is_copyable (a)
  return rawget(getmetatable(a) or {}, '__copy')
end

local function is_callable (a)
  return type(a) == 'function' or rawget(getmetatable(a) or {}, '__call')
end

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

local function var2obj (self) -- convert variables to an object
  local sv = self[var]
  sv[var] = {}                                -- set var
  sv.name = rawget(self,name)                 -- set name
  return setmetatable(sv, getmetatable(self)) -- set parent
end

local function proxy (self, k, v) -- wrap table with a proxy
  local sv = self[var]
  sv[k] = setmetatable({name=k, [var]=v, [obj]=rawget(self,obj) or self}, MT)
  return sv[k]
end

local function eval_tbl(self, k, v) -- read-eval table members
  return is_function(v) and v(self[obj]) or
         is_table   (v) and proxy(self, k, v) or v
end

local function eval_obj(self, k, v) -- read-eval object members
  return is_function(v) and v(self) or
         is_table   (v) and proxy(self, k, v) or v
end

-- proxy for controlling tables stored in object variables

function MT:__index (k) -- (+__eval)
  local v = self[var][k]
--  io.write('++ ', rawget(self[obj], 'name') or tostring(self[obj]), '.', self.name, '[', tostring(k), '] = ', tostring(v), '\n')
  return v and eval_tbl(self, k, v)
end

function MT:__newindex (k, nv) -- (+__copy)
  local v = self[var][k]
  local c = is_copyable(v) or is_copyable(nv)
  self[var][k] = c and c(v,nv) or nv
end

function MT:__len    () return #self[var]        end
function MT:__pairs  () return pairs(self[var])  end
function MT:__ipairs () return ipairs(self[var]) end

-- proxy for controlling object variables

function M:__call (a) -- object ctor
  if is_table(a) then
    if self[var] == var then self[var] = a; return self end
    return setmetatable( {        [var]=a  }, obj2cls(self) )
  elseif is_string(a) then
    return setmetatable( {name=a, [var]=var}, obj2cls(self) )
  end
  error("invalid object argument")
end

function M:__index (k) -- (+__eval+inheritance)
  local v = self[var][k]
--  io.write('++ ', rawget(self, 'name') or tostring(self), '.', tostring(k), " = ", tostring(v), "\n")
  return v and eval_obj(self, k, v) or getmetatable(self)[k]
end

function M:__newindex (k, nv) -- (+__copy)
  local v = self[var][k]
  local c = is_copyable(v) or is_copyable(nv)
  self[var][k] = c and c(v,nv) or nv
end

function M:__len     () return #self[var]         end
function M:__pairs   () return pairs(self[var])   end
function M:__ipairs  () return ipairs(self[var])  end
function M:parent    () return getmetatable(self) end
function M:variables () return self[var]          end

function M:get (k) -- idem __index (-__eval)
  return self[var][k] or getmetatable(self)[k]
end

function M:set (a, v) -- idem __newindex (+shallow copy)
  local sv = self[var]
  if v ~= nil then
    sv[a] = v
  else -- shallow copy (slow)
    for k,v in pairs(a) do sv[k] = v end
  end
  return self
end

function M:set_function (k, f)
  assert(is_callable(f), "invalid set_function argument")
  rawset(self,k,f)
  return self
end

function M:make_class ()
  assert(self[var] ~= nil, "invalid or incomplete object")
  return self:var2obj():init_class()
end

-- debug

function M:dump (file, level, indent)
  local fp = file or io.stdout
  local lv = level or 1
  local id = indent or 1
  local sv = self[var]
  local pa = self:parent()
  if id == 1 then -- header
    fp:write("objdump '", self.name, "' [", tostring(sv), "]\n")
  end
  for k,v in pairs(sv) do -- keys
    for i=1,id do fp:write("  ") end -- indent
    fp:write(tostring(k), ": ",
             tostring(getmetatable(v) == MT and v[var] or v), "\n")
  end
  if lv > 1 and pa ~= M then -- parent
    for i=1,id do fp:write("  ") end -- indent
    fp:write("parent '", pa.name, "' [", tostring(pa[var]), "]\n")
    pa:dump(fp, lv-1, id+1)
  end
  fp:write("\n")
end

------------------------------------------------------------------------------o
return setmetatable( {name='Object', [var]={}}, M ) -- root
