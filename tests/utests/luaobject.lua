--[=[
 o-----------------------------------------------------------------------------o
 |
 | Object model (pure Lua) regression tests
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - Provide regression test suites for the object model without extension.

 o-----------------------------------------------------------------------------o
]=]

-- expected from other modules ------------------------------------------------o

require 'table.clear'
require 'table.new'

-- operators
local bit = require 'bit'
local tobit = bit.tobit
local bnot, band, bor, rol = bit.bnot, bit.band, bit.bor, bit.rol
local lshift, rshift  = bit.lshift, bit.rshift

local bset = function (x,n) return bor (x, rol( 1, n))      end
local btst = function (x,n) return band(x, rol( 1, n)) ~= 0 end
local bclr = function (x,n) return band(x, rol(-2, n))      end
local lt   = function (x,y) return x <  y                   end
local le   = function (x,y) return x <= y                   end
local gt   = function (x,y) return x >  y                   end
local ge   = function (x,y) return x >= y                   end

-- metatables & metamethods
local ffi = require 'ffi'
local typeof, miscmap in ffi
assert(miscmap, "missing MAD extension (no cdata metatable access)")

local metaname = {
  '__add', '__call', '__concat', '__copy', '__div', '__eq', '__exec', '__gc',
  '__index', '__init', '__ipairs', '__le', '__len', '__lt', '__metatable',
  '__mod', '__mode', '__mul', '__new', '__newindex', '__pairs', '__pow',
  '__same', '__sub', '__tostring', '__totable', '__unm',
}
for _,v in ipairs(metaname) do metaname[v] = v end -- build dict of metaname

local function invalid_use (a)
  error("invalid use of object <" .. tostring(a) .. ">", 2)
end

local function get_metatable (a)
  return type(a) == 'cdata' and (a.__metatable or miscmap[-tonumber(typeof(a))])
         or getmetatable(a)
end

local function has_metatable (a)
  return not not get_metatable(a)
end

local function get_metamethod (a, f)
  local mt = get_metatable(a)
  return mt and rawget(mt,f)
end

local function has_metamethod (a, f)
  local mt = get_metatable(a)
  local mm = mt and rawget(mt,f)
  return not not mm and mm ~= invalid_use
end

local function has_metamethod_ (a, f)
  local mt = get_metatable(a)
  local mm = mt and rawget(mt,f)
  return not not mm -- false or nil -> false
end

-- types

local function is_nil (a)
  return type(a) == 'nil'
end

local function is_boolean (a)
  return type(a) == 'boolean'
end

local function is_number (a)
  return type(a) == 'number'
end

local function is_string (a)
  return type(a) == 'string'
end

local function is_function (a)
  return type(a) == 'function'
end

local function is_table (a)
  return type(a) == 'table'
end

local function is_rawtable (a)
  return type(a) == 'table' and rawequal(getmetatable(a),nil)
end

local function is_metaname (a)
  return metaname[a] == a
end

-- concepts

local function is_callable (a)
  return type(a) == 'function' or has_metamethod_(a,'__call')
end

local function is_iterable (a)
 return type(a) == 'table' or has_metamethod (a,'__ipairs')
end

local function is_mappable (a)
 return type(a) == 'table' or has_metamethod (a,'__pairs')
end

-- iterator: pairs = ipairs + kpairs

local function kpairs_iter (tbl, key)
  local k, v = key
  repeat k, v = tbl.nxt(tbl.dat, k)                -- discard ipairs indexes
  until type(k) ~= 'number' or (k%1) ~= 0 or k > tbl.lst or k < 1
  return k, v
end

local function kpairs (tbl, lst_)
  assert(is_mappable(tbl), "invalid argument #1 (mappable expected)")
  local nxt, dat, ini = pairs(tbl)
  local lst = lst_
  if not lst and is_iterable(tbl) then
    for i in ipairs(tbl) do lst = i end
  end
  if not lst then return nxt, dat, ini end
  return kpairs_iter, { nxt=nxt, dat=dat, lst=lst }, ini
end

-- searching

local function bsearch (tbl, val, cmp_, low_, high_)
  assert(is_iterable(tbl), "invalid argument #1 (iterable expected)")
  if is_number(cmp_) and is_nil(high_) then
    cmp_, low_, high_ = nil, cmp_, low_ -- right shift
  end
  if not (is_nil(cmp_) or is_callable(cmp_)) then
    error("invalid argument #3 (callable expected)")
  end
  local cmp  = cmp_ or lt -- cmp must be the same used by table.sort
  local low  = is_number(low_ ) and max(low_ , 1   ) or 1
  local high = is_number(high_) and min(high_, #tbl) or #tbl
  local cnt, mid, stp, tst = high-low+1
  while cnt > 0 do
    stp = rshift(cnt,1)
    mid = low+stp
    tst = cmp(tbl[mid], val)
    low = tst and mid+1 or low
    cnt = tst and rshift(cnt-1,1) or stp
  end
  return low -- tbl[low] <= val
end

-- functors

local _fun, _fun2 = {}, {}
local fct_mt, fct_mtc

local function functor (f)
  assert(is_callable(f), "invalid argument #1 (callable expected)")
  return setmetatable({[_fun]=f}, fct_mt)
end

local function compose (f, g)
  assert(is_callable(f), "invalid argument #1 (callable expected)")
  assert(is_callable(g), "invalid argument #2 (callable expected)")
  return setmetatable({[_fun]=f, [_fun2]=g}, fct_mtc)
end

local function is_functor (a)
  return is_table(a) and rawget(a,_fun) ~= nil
end

local str = function(s) return string.format("functor: %p", s) end
local err = function()  error("forbidden access to functor", 2) end

fct_mt = {
  __pow       = compose,
  __call      = function(s,...) return rawget(s,_fun)(...) end,
  __index     = function(s,k)   return s(k)                end,
  __tostring  = str,
  __len = err, __newindex = err, __ipairs = err, __pairs = err,
}

fct_mtc = {
  __pow       = compose,
  __call      = function(s,...) return rawget(s,_fun)(rawget(s,_fun2)(...)) end,
  __index     = function(s,k)   return s(k)                                 end,
  __tostring  = str,
  __len = err, __newindex = err, __ipairs = err, __pairs = err,
}

-- implementation -------------------------------------------------------------o

local M = {}

-- Root of all objects, forward declaration
local object

-- object members
local _var = {} -- hidden key

-- reserved flags (bits)
local _flg = {} -- hidden key
local flg_ro, flg_cl = 0, 1 -- flags id for readonly and class
local flg_free = 2          -- used flags (0 - 1), free flags (2 - 31)

-- instance and metatable of 'incomplete objects' proxy
local var0 = setmetatable({}, {
  __index     = function() error("forbidden read access to incomplete object" , 2) end,
  __newindex  = function() error("forbidden write access to incomplete object", 2) end,
  __metatable = false,
})

-- helpers

local function name (a)
  local var = rawget(a,_var)
  return rawget(var,'__id') or ('? <: ' .. var.__id)
end

local function init (a)
  local init = rawget(getmetatable(a), '__init')
  if init then return init(a) end
  return a
end

local function parent (a)
  return getmetatable(rawget(a,'__index'))
end

local function fclass (a)
  return btst(rawget(a,_flg) or 0, flg_cl)
end

local function freadonly (a)
  return btst(rawget(a,_flg) or 0, flg_ro)
end

local function set_class (a)
  rawset(a,_flg, bset(rawget(a,_flg) or 0, flg_cl))
  return a
end

local function is_object (a) -- exported
  return is_table(a) and rawget(a,_var) ~= nil
end

local function is_class (a) -- exported
  return is_table(a) and fclass(a)
end

local function is_instanceOf (a, b) -- exported
  if is_object(a) and is_class(b) then
    repeat a = parent(a) until not a or rawequal(a,b)
    return not not a
  end
  return false
end

-- metamethods

local MT = {}

-- objects are proxies controlling variables access and inheritance
function MT:__call (a, b) -- object constructor (define the object-model)
  if is_string(a) or is_nil(a) then                     -- [un]named object
    if is_nil(b) then
      local obj = {__id=a, [_var]=var0, __index=rawget(self,_var)} -- proxy
      return setmetatable(obj, getmetatable(self))      -- incomplete object
    elseif is_rawtable(b) then
      local obj = {[_var]=b, __index=rawget(self,_var)} -- proxy
      b.__id=a ; setmetatable(b, obj) ; set_class(self) -- set fast inheritance
      return init(setmetatable(obj, getmetatable(self)))-- complete object
    end
  elseif is_rawtable(a) then
    if rawget(self,_var) == var0 then                   -- finalize named object
      a.__id, self.__id = self.__id, nil
      rawset(self,_var, setmetatable(a, self));         -- set fast inheritance
      set_class(parent(self))
      return init(self)
    else                                                -- unnamed object
      local obj = {[_var]=a, __index=rawget(self,_var)} -- proxy
      setmetatable(a, obj) ; set_class(self)            -- set fast inheritance
      return init(setmetatable(obj, getmetatable(self)))-- complete object
    end
  end
  error(is_nil(b) and "invalid argument #1 (string or raw table expected)"
                  or  "invalid argument #2 (raw table expected)", 2)
end

local function raw_len (self)
  return rawlen(rawget(self,_var))    -- no inheritance
end

local function raw_get (self, k)
  return rawget(rawget(self,_var),k)  -- no inheritance nor function evaluation
end

local function raw_set (self, k, v)
  rawset(rawget(self,_var), k, v)     -- no protection!!
end

local function var_raw (self, k)
  return rawget(self,_var)[k]         -- no function evaluation with inheritance
end

local function var_val (self, k, v)   -- string key with value function
  if type(k) == 'string' and type(v) == 'function'
  then return v(self)
  else return v end
end

local function var_get (self, k) -- reusing var_raw and var_val kills the inlining
  local v = rawget(self,_var)[k]
  if type(k) == 'string' and type(v) == 'function'
  then return v(self)
  else return v end
end

function MT:__index (k)          -- reusing var_raw and var_val kills the inlining
  local v = rawget(self,_var)[k]
  if type(k) == 'string' and type(v) == 'function'
  then return v(self)
  else return v end
end

function MT:__newindex (k, v)
  if freadonly(self) or type(k) == 'string' and string.sub(k,1,2) == '__' then
    error("forbidden write access to '" .. name(self) ..
          "' (readonly object or variable)", 2)
  end
  rawget(self,_var)[k] = v      -- note: must use [k] for var0
end

function MT:__len ()
  local var = rawget(self,_var)
  if is_nil(var[1]) then return 0 end -- fast
  while is_nil(rawget(var,1)) do      -- slow
    var  = rawget(self,'__index')
    self = getmetatable(var)
  end
  return rawlen(var)
end

local function iter (var, key) -- scan only numbers and strings
  local k, v = next(var, key)
  while type(k) ~= 'string' and type(k) ~= 'number' and k do
    k, v = next(var, k)
  end
  return k, v
end

local function pairs_iter (self)
  return iter, rawget(self,_var), nil
end

local function ipairs_iter (self)
  return ipairs(rawget(self,_var))
end

MT.__pairs  =  pairs_iter
MT.__ipairs = ipairs_iter

function MT:__tostring()
  return string.format("object: '%s' %p", name(self), self)
end

-- methods

local function is_readonly (self)
  assert(is_object(self), "invalid argument #1 (object expected)")
  return freadonly(self)
end

local function set_readonly (self, set_)
  assert(is_object(self), "invalid argument #1 (object expected)")
  if set_ ~= false
  then rawset(self, _flg, bset(rawget(self,_flg) or 0, flg_ro))
  else rawset(self, _flg, bclr(rawget(self,_flg) or 0, flg_ro))
  end
  return self
end

local function get_variables (self, lst, noeval_)
  assert(is_object(self) , "invalid argument #1 (object expected)")
  assert(is_iterable(lst), "invalid argument #2 (iterable expected)")
  local n   = #lst
  local res = table.new(0,n)
  local get = noeval_ == true and var_raw or var_get
  for i=1,n do res[lst[i]] = get(self, lst[i]) end
  return res -- key -> val
end

local function set_variables (self, tbl, override_)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(is_mappable(tbl)   , "invalid argument #2 (mappable expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  local id  = rawget(var,'__id')
  for k,v in pairs(tbl) do
    assert(is_nil(rawget(var,k)) or override_ ~= false, "cannot override variable")
    var[k] = v
  end
  var.__id = id
  return self
end

local function wrap_variables (self, tbl)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(is_mappable(tbl)   , "invalid argument #2 (mappable expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  local id  = rawget(var,'__id')
  for k,f in pairs(tbl) do
    local v, newv = var[k]
    assert(not is_nil(v) , "invalid variable (nil value)")
    assert(is_callable(f), "invalid wrapper (callable expected)")
    if is_callable(v) then newv = f(v) else newv = f(\ v) end -- simplify user's side.
    if is_functor(v) and not is_functor(newv) then
      newv = functor(newv)                   -- newv must maintain v's semantic.
    end
    var[k] = newv
  end
  var.__id = id
  return self
end

local function set_methods (self, tbl, override_, strict_)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(is_mappable(tbl)   , "invalid argument #2 (mappable expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  local id  = rawget(var,'__id')
  for k,f in pairs(tbl) do
    assert(is_string(k), "invalid key (function name expected)")
    assert(is_callable(f) or strict_ == false, "invalid value (callable expected)")
    assert(is_nil(rawget(var,k)) or override_ ~= false, "cannot override function")
    var[k] = is_function(f) and functor(f) or f
  end
  var.__id = id
  return self
end

local function set_metamethods (self, tbl, override_, strict_)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(is_mappable(tbl)   , "invalid argument #2 (mappable expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local sm, pm = getmetatable(self), getmetatable(parent(self)) or MT
  if sm == pm then -- create new metatable if same as parent
    assert(not fclass(self), "invalid metatable write access (unexpected class)")
    sm=table.new(0,8) for k,v in pairs(pm) do sm[k] = v end
    pm.__metatable = nil -- unprotect change
    setmetatable(self, sm)
    pm.__metatable, sm.__metatable = pm, sm
  end
  for k,mm in pairs(tbl) do
    assert(is_metaname(k) or strict_ == false, "invalid key (metamethod expected)")
    assert(is_nil(rawget(sm,k)) or override_ == true, "cannot override metamethod")
    sm[k] = mm
  end
  return self
end

local function final_err (self)
  error("invalid object creation ('"..name(self).."' is qualified as final)", 2)
end

local function set_final (self)
  return set_metamethods(self, {__call=final_err}, true)
end

local function get_varkeys (self, class_)
  assert(is_object(self)                    , "invalid argument #1 (object expected)")
  assert(is_nil(class_) or is_object(class_), "invalid argument #2 (object expected)")
  local lst, key = table.new(8,1), table.new(0,8)
  while self and not rawequal(self, class_) do
    for k,v in kpairs(self) do
      if not (key[k] or is_functor(v)) and is_string(k) and string.sub(k,1,2) ~= '__'
      then lst[#lst+1], key[k] = k, k
      end
    end
    self = parent(self)
  end
  assert(rawequal(self, class_),"invalid argument #2 (parent of argument #1 expected)")
  return lst
end

local function insert (self, idx_, val)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  table.insert(rawget(self,_var), idx_, val)
  return self
end

local function remove (self, idx_)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  return table.remove(rawget(self,_var), idx_)
end

local function move (self, idx1, idx2, idxto, dest_)
  dest_ = dest_ or self
  assert(is_object(self)     , "invalid argument #1 (object expected)")
  assert(is_object(dest_)    , "invalid argument #2 (object expected)")
  assert(not freadonly(dest_), "forbidden write access to readonly object")
  table.move(rawget(self,_var), idx1, idx2, idxto, rawget(dest_,_var))
  return dest_
end

local function sort (self, cmp_)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  table.sort(rawget(self,_var), cmp_)
  return self
end

local function bsearch_ (self, val, cmp_, low_, high_)
  assert(is_object(self), "invalid argument #1 (object expected)")
  return bsearch(rawget(self,_var), val, cmp_, low_, high_)
end

local function lsearch_ (self, val, equ_, low_, high_)
  assert(is_object(self), "invalid argument #1 (object expected)")
  return lsearch(rawget(self,_var), val, equ_, low_, high_)
end

local function clear_array (self)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  for i=1,rawlen(var) do var[i]=nil end
  return self
end

local function clear_variables (self)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  local id  = rawget(var,'__id')
  for k in kpairs(self) do var[k]=nil end
  var.__id = id
  return self
end

local function clear_all (self)
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(not freadonly(self), "forbidden write access to readonly object")
  local var = rawget(self,_var)
  local id  = rawget(var,'__id')
  for k in pairs_iter(self) do var[k]=nil end -- table.clear destroys all keys
  var.__id = id
  return self
end

-- inheritance

local function set_parent (self, newp)
  assert(is_object(self), "invalid argument #1 (object expected)")
  assert(is_object(newp), "invalid argument #2 (object expected)")
  if freadonly(self) then
    error("forbidden write access to readonly object '" .. name(self) .. "'", 2)
  end
  local spar = self.parent
  if getmetatable(newp) ~= getmetatable(spar) then
    error("new and current parent do not share same metamethods")
  end
  if newp.parent ~= spar.parent then
    error("new and current parent do not inherit from same direct parent")
  end
  rawset(self,'__index', rawget(newp,_var))
  set_class(newp)
  return self
end

-- copy

local function same (self, name_)
  assert(is_object(self)                  ,"invalid argument #1 (object expected)")
  assert(is_nil(name_) or is_string(name_),"invalid argument #2 (string expected)")
  -- same shares the same parent
  local par = parent(self)
  local sam = par(name_, {})
  -- metatable
  local sm, pm = getmetatable(self), getmetatable(par)
  if sm ~= pm then -- copy metatable
    local cm=table.new(0,8) for k,v in pairs(sm) do cm[k] = v end
    sm.__metatable = nil
    setmetatable(sam, cm)
    sm.__metatable, cm.__metatable = sm, cm
  end
  return sam
end

local function copy (self, name_)
  assert(is_object(self)                  ,"invalid argument #1 (object expected)")
  assert(is_nil(name_) or is_string(name_),"invalid argument #2 (string expected)")
  local cpy = same(self, name_ or raw_get(self,'__id'))
  local var, cvar = rawget(self,_var), rawget(cpy,_var)
  local id  = rawget(cvar,'__id')
  for k,v in pairs_iter(self) do cvar[k] = v end
  cvar.__id = id
  return cpy
end

MT.__same = same
MT.__copy = copy

-- flags

local flg_mask = lshift(-1, flg_free)
local flg_notmask = bnot(flg_mask)

local function test_flag (self, n)
  assert(is_object(self), "invalid argument #1 (object expected)")
  assert(is_number(n)   , "invalid argument #2 (number expected)")
  return btst(rawget(self,_flg) or 0, n)
end

local function set_flag (self, n)
  assert(is_object(self), "invalid argument #1 (object expected)")
  assert(is_number(n)   , "invalid argument #2 (number expected)")
  if n >= flg_free then
    rawset(self, _flg, bset(rawget(self,_flg) or 0, n))
  end
  return self
end

local function clear_flag (self, n)
  assert(is_object(self), "invalid argument #1 (object expected)")
  assert(is_number(n)   , "invalid argument #2 (number expected)")
  if n >= flg_free then
    rawset(self, _flg, bclr(rawget(self,_flg) or 0, n))
  end
  return self
end

local function get_flags (self)
  assert(is_object(self), "invalid argument #1 (object expected)")
  return rawget(self,_flg) or 0
end

local function test_flags (self, flags)
  assert(is_object(self) , "invalid argument #1 (object expected)")
  assert(is_number(flags), "invalid argument #2 (number expected)")
  return band(rawget(self,_flg) or 0, flags) ~= 0
end

local function set_flags (self, flags)
  assert(is_object(self) , "invalid argument #1 (object expected)")
  assert(is_number(flags), "invalid argument #2 (number expected)")
  flags = band(flags, flg_mask)
  rawset(self, _flg, bor(rawget(self, _flg) or 0, flags))
  return self
end

local function clear_flags (self, flags)
  assert(is_object(self) , "invalid argument #1 (object expected)")
  assert(is_number(flags), "invalid argument #2 (number expected)")
  flags = band(flags, flg_mask)
  rawset(self, _flg, band(rawget(self, _flg) or 0, bnot(flags)))
  return self
end

-- environments

local _env = {} -- hidden key

local function open_env (self, ctx_)
  assert(is_object(self), "invalid argument #1 (object expected)")
  assert(is_nil(ctx_) or is_function(ctx_) or is_number(ctx_) and ctx_ >= 1,
                          "invalid argument #2 (not a function or < 1)")
  ctx_ = is_function(ctx_) and ctx_ or is_number(ctx_) and ctx_+1 or 2
  assert(is_nil(rawget(self,_env)), "invalid environment (already open)")
  rawset(self, _env, { ctx=ctx_, env=getfenv(ctx_) })
  rawset(self, self.__id, self) -- self reference
  setfenv(ctx_, self)
  return self
end

local function is_open_env (self)
  assert(is_object(self), "invalid argument #1 (object expected)")
  return not is_nil(rawget(self,_env))
end

local function reset_env (self) -- if an error occurs while in the environment
  assert(is_object(self), "invalid argument #1 (object expected)")
  rawset(self, _env, nil)
  rawset(self, self.__id, nil) -- clear self reference
  return self
end

local function close_env (self)
  assert(is_object(self), "invalid argument #1 (object expected)")
  local env = rawget(self,_env)
  assert(not is_nil(env), "invalid environment (not open)")
  setfenv(env.ctx, env.env)
  return reset_env(self)
end

local function dump_env (self) -- for debug
  for k,v in pairs(rawget(self,_var)) do
    if is_rawtable(v) then
      for k,v in pairs(v) do
        print(k,'=',v)
      end
    elseif is_object(v) then
      print(k,'=',name(v))
    else
      print(k,'=',v)
    end
  end
end

-- I/O ------------------------------------------------------------------------o

-- dump obj members (including controlled inheritance)
local function dumpobj (self, filnam_, class_, pattern_)
  if is_object(filnam_) and is_nil(pattern_) then
    filnam_, class_, pattern_ = nil, filnam_, class_ -- right shift
  end
  if is_string(class_) and is_nil(pattern_) then
    class_, pattern_ = nil, class_                   -- right shift
  end

  class_, pattern_ = class_ or object, pattern_ or ''
  assert(is_object(self)    , "invalid argument #1 (object expected)")
  assert(is_object(class_)  , "invalid argument #3 (object expected)")
  assert(is_string(pattern_), "invalid argument #4 (string expected)")

  local tostring in MAD
  local n, cnt, res, spc, str = 0, {}, {}, ""
  while self and not rawequal(self, class_) do
    local var = rawget(self,_var)
    -- header
    local id = rawget(var, '__id')
    n, str = n+1, id and (" '" .. id .. "'") or ""
    res[n] = spc .. "+ " .. tostring(self)
    spc = spc .. "   "
    -- variables
    for k,v in kpairs(self) do
      if is_string(k) and string.sub(k,1,2) ~= '__' and string.find(k,pattern_) then
        str = spc .. tostring(k)
        if is_string(v) then
          str = str .. " : '" .. tostring(v):sub(1,15) .. "'"
        elseif is_function(v) then
          str = str .. " := " .. tostring(v(self))
        else
          str = str .. " :  " .. tostring(v)
        end
        if cnt[k]
        then str = str .. " (" .. string.rep('*', cnt[k]) .. ")" -- mark overrides
        else cnt[k] = 0
        end
        cnt[k], n = cnt[k]+1, n+1
        res[n] = str
      end
    end
    self = parent(self)
  end
  assert(rawequal(self, class_), "invalid argument #2 (parent of argument #1 expected)")

  -- return result as a string
  if filnam_ == '-' then return table.concat(res, '\n') end

  -- dump to file
  local file = openfile(filnam_, 'w', '.dat')
  for _,s in ipairs(res) do file:write(s,'\n') end
  if is_string(filnam_) then file:close() else file:flush() end

  return self
end

-- members --------------------------------------------------------------------o

M.__id  = 'object'
M.__par = parent
M.first_free_flag = flg_free

-- methods
M.is_class        = functor( is_class        )
M.is_readonly     = functor( is_readonly     )
M.is_instanceOf   = functor( is_instanceOf   )

M.set_parent      = functor( set_parent      )
M.set_readonly    = functor( set_readonly    )
M.set_final       = functor( set_final       )

M.get_varkeys     = functor( get_varkeys     )
M.get_variables   = functor( get_variables   )
M.set_variables   = functor( set_variables   )
M.wrap_variables  = functor( wrap_variables  )

M.set_methods     = functor( set_methods     )
M.set_metamethods = functor( set_metamethods )

M.insert          = functor( insert          )
M.remove          = functor( remove          )
M.move            = functor( move            )
M.sort            = functor( sort            )
M.bsearch         = functor( bsearch_        )
M.lsearch         = functor( lsearch_        )
M.clear_array     = functor( clear_array     )
M.clear_variables = functor( clear_variables )
M.clear_all       = functor( clear_all       )

M.same            = functor( same            )
M.copy            = functor( copy            )

M.raw_len         = functor( raw_len         )
M.raw_get         = functor( raw_get         )
M.raw_set         = functor( raw_set         )

M.var_raw         = functor( var_raw         )
M.var_val         = functor( var_val         )
M.var_get         = functor( var_get         )

M.set_flag        = functor( set_flag        )
M.test_flag       = functor( test_flag       )
M.clear_flag      = functor( clear_flag      )
M.get_flags       = functor( get_flags       )
M.set_flags       = functor( set_flags       )
M.test_flags      = functor( test_flags      )
M.clear_flags     = functor( clear_flags     )

M.open_env        = functor( open_env        )
M.reset_env       = functor( reset_env       )
M.close_env       = functor( close_env       )
M.is_open_env     = functor( is_open_env     )

M.dumpobj         = functor( dumpobj         )

-- aliases
M.parent = parent
M.name   = function(s) return s.__id end
M.set    = M.set_variables
M.get    = M.get_variables

-- metatables -----------------------------------------------------------------o

-- root object variables = module
object = setmetatable({[_var]=M}, MT)

 -- parent link
setmetatable(M, object)

-- protect against changing metatable
MT.__metatable = MT

-- set as readonly
object:set_readonly()

-- end of object model --------------------------------------------------------o

local utest
if MAD == nil then
  utest = require("luaunit")
else
  utest = MAD.utest
end

local assertEquals = utest.assertEquals
local assertAlmostEquals = utest.assertAlmostEquals
local assertErrorMsgContains = utest.assertErrorMsgContains
local assertNil = utest.assertNil
local assertTrue = utest.assertTrue
local assertFalse = utest.assertFalse
local assertNotEquals = utest.assertNotEquals
local assertNotNil = utest.assertNotNil
local assertStrContains = utest.assertStrContains

-- regression test suite ------------------------------------------------------o

TestLuaObject = {}
TestLuaObjectErr = {}
local myFunc = function() end
local objectErr = {1, true, '', {}, myFunc}
local _msg = {
  "invalid argument #1 (object expected)",
  "invalid argument #2 (object expected)",
  "forbidden write access to readonly object"
}

function TestLuaObject:testConstructor()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object
  local p3 = object('p3',{})

  local p00 = p0 'p00' {}
  local p01 = p0 {}
  local p03 = p0
  local p04 = p0('p04',{})

  -- ctor equivalence
  assertEquals(p04, object 'p04' {})
  assertEquals(p04, p0 'p04' {})

  -- read
  assertNil( p0.a )
  assertNil( p1.a )
  assertNil( p2.a )
  assertNil( p3.a )

  -- write
  p0.a = ''   assertEquals( p0.a, '' )
  p1.a = ''   assertEquals( p1.a, '' )
  p3.a = ''   assertEquals( p3.a, '' )

  -- read child
  assertEquals( p00.a, '' )
  assertEquals( p01.a, '' )
  assertEquals( p03.a, '' )
  assertEquals( p04.a, '' )

  -- write child
  p00.a = '0'   assertEquals( p00.a, '0' )
  p01.a = '0'   assertEquals( p01.a, '0' )
  p03.a = '0'   assertEquals( p03.a, '0' )
  p04.a = '0'   assertEquals( p04.a, '0' )
end

function TestLuaObjectErr:testConstructor()
  local p0 = object 'p0' {}
  local p2 = object 'p2'
  local p3 = object
  local p00 = p0 'p02'

  local notraw_table = setmetatable({}, {})

  local get = function(s,k) return s[k] end
  local set = function(s,k,v) s[k]=v end
  local msg = {
    "forbidden read access to incomplete object",
    "forbidden write access to incomplete object",
    "invalid argument #1 (string or raw table expected)",
    "invalid argument #2 (raw table expected)",
    "forbidden write access to 'object' (readonly object or variable)",
  }
  local a
  assertErrorMsgContains(msg[1], get, p2, a)      -- read
  assertErrorMsgContains(msg[1], get, p00, a)     -- read child

  assertErrorMsgContains(msg[2], set, p2, a, '')  -- write
  assertErrorMsgContains(msg[2], set, p00, a, '') -- write child

  assertErrorMsgContains(msg[5], set, p3, a, '')  -- write

  assertErrorMsgContains(msg[3], object, true)
  assertErrorMsgContains(msg[3], object, 1)
  assertErrorMsgContains(msg[3], object, myFunc)
  assertErrorMsgContains(msg[3], object, object)
  assertErrorMsgContains(msg[3], object, notraw_table)

  assertErrorMsgContains(msg[4], object, 'p', true)
  assertErrorMsgContains(msg[4], object, 'p', 1)
  assertErrorMsgContains(msg[4], object, 'p', '1')
  assertErrorMsgContains(msg[4], object, 'p', myFunc)
  assertErrorMsgContains(msg[4], object, 'p', object)
  assertErrorMsgContains(msg[4], object, 'p', notraw_table)
end

function TestLuaObject:testInheritance()
  local p0 = object {}
  local p1 = p0 { x=3, y=2, z=1  }
  local p2 = p1 { x=2, y=1 }
  local p3 = p2 { x=1  }
  local p4 = p3 { }
  local vs = {'x','y','z'}

  assertEquals   ( p0:get(vs), {} )
  assertEquals   ( p0        , {} )

  assertEquals   ( p1:get(vs), { x=3, y=2, z=1 } )
  assertEquals   ( p1        , { x=3, y=2, z=1 } )
  assertNotEquals( p1        , { x=3, y=2 } )

  assertEquals   ( p2:get(vs), { x=2, y=1, z=1 } )
  assertEquals   ( p2        , { x=2, y=1 } )
  assertNotEquals( p2        , { x=2 } )

  assertEquals   ( p3:get(vs), { x=1, y=1, z=1 } )
  assertEquals   ( p3        , { x=1 } )
  assertNotEquals( p3        , { x=1, y=1 } )

  assertEquals   ( p4:get(vs), { x=1, y=1, z=1 } )
  assertEquals   ( p4        , {} )
  assertNotEquals( p4        , { x=1 } )

  assertEquals   ( {p1.x, p1.y, p1.z}, {3,2,1})
  assertEquals   ( {p2.x, p2.y, p2.z}, {2,1,1})
  assertEquals   ( {p3.x, p3.y, p3.z}, {1,1,1})
  assertEquals   ( {p4.x, p4.y, p4.z}, {1,1,1})

  p2:set{x=5, y=6}  p4:set{y=5, z=6}

  assertEquals   ( p0:get(vs), {} )
  assertEquals   ( p0        , {} )

  assertEquals   ( p1:get(vs), { x=3, y=2, z=1 } )
  assertEquals   ( p1        , { x=3, y=2, z=1 } )
  assertNotEquals( p1        , { x=3, y=2 } )

  assertEquals   ( p2:get(vs), { x=5, y=6, z=1 } )
  assertEquals   ( p2        , { x=5, y=6 } )
  assertNotEquals( p2        , { x=5 } )
  assertNotEquals( p2        , { x=2, y=1 } )

  assertEquals   ( p3:get(vs), { x=1, y=6, z=1 } )
  assertEquals   ( p3        , { x=1 } )
  assertNotEquals( p3        , { x=1, y=6 } )
  assertNotEquals( p3        , { x=1, y=1 } )

  assertEquals   ( p4:get(vs), { x=1, y=5, z=6 } )
  assertEquals   ( p4        , { y=5, z=6 } )
  assertNotEquals( p4        , { x=1 } )

  assertEquals   ( {p1.x, p1.y, p1.z}, {3,2,1})
  assertEquals   ( {p2.x, p2.y, p2.z}, {5,6,1})
  assertEquals   ( {p3.x, p3.y, p3.z}, {1,6,1})
  assertEquals   ( {p4.x, p4.y, p4.z}, {1,5,6})
end

function TestLuaObject:testIsObject()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object 'p2'
  local p3 = object

  local p00 = p0 'p00' {}
  local p01 = p0 {}
  local p02 = p0 'p02'
  local p03 = p0

  assertTrue ( is_object(p0) )
  assertTrue ( is_object(p1) )
  assertTrue ( is_object(p2) )
  assertTrue ( is_object(p3) )
  assertTrue ( is_object(p00) )
  assertTrue ( is_object(p01) )
  assertTrue ( is_object(p02) )
  assertTrue ( is_object(p03) )
  assertFalse( is_object(nil)  )
  assertFalse( is_object(1)    )
  assertFalse( is_object({})   )
  assertFalse( is_object("yes"))
  assertFalse( is_object(myFunc))
end

function TestLuaObject:testIsClass()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object 'p2'
  local p3 = object
  local p00 = p0 'p00' {}
  local p01 = p0 {}
  local p02 = p0 'p02'
  local p03 = p0

  assertTrue ( is_class(object) )
  assertTrue ( is_class(p0) )
  assertTrue ( is_class(p3) )
  assertTrue ( is_class(p03) )
  assertFalse( is_class(p00) )
  assertFalse( is_class(p01) )
  assertFalse( is_class(p02) )
  assertFalse( is_class(p1)  )
  assertFalse( is_class(p2)  )
  assertFalse( is_class(nil)  )
  assertFalse( is_class(1)    )
  assertFalse( is_class({})   )
  assertFalse( is_class("yes"))
  assertFalse( is_class(myFunc))
end

function TestLuaObject:testIsReadonly()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object 'p2' {}
  local p3 = object
  local p00 = p0 'p00' {}
  local p01 = p0 {}
  local p02 = p0 'p02' {}
  local p03 = p0

  assertFalse( p0 :is_readonly() )
  assertFalse( p1 :is_readonly() )
  assertFalse( p2 :is_readonly() )
  assertTrue ( p3 :is_readonly() )
  assertFalse( p00:is_readonly() )
  assertFalse( p01:is_readonly() )
  assertFalse( p02:is_readonly() )
  assertFalse( p03:is_readonly() )
  assertTrue ( p03:set_readonly(true):is_readonly() )
  p0:set_readonly(true)
  assertTrue ( p0:is_readonly()  )
  assertFalse( p00:is_readonly() )
end

function TestLuaObject:testSetParent()
  local p0 = object 'p0' {}
  local p00 = p0  "p00" {}
  local p01 = p0  "p01" { a = true }
  local p02 = p00 "p02" { b = true }

  assertNil(p02.a)
  p02:set_parent(p01)
  assertTrue(p02.a)
  assertTrue(p02.b)
end

function TestLuaObjectErr:testSetParent()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }:set_readonly()

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set_parent, objectErr[i])
    assertErrorMsgContains(_msg[2], p0.set_parent, p0, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.set_parent, p1, p1)
  assertErrorMsgContains(_msg[3], object.set_parent, object, p1)
end

function TestLuaObject:testIsInstanceOf()
  local p0 = object {}
  local p1 = p0 { }
  local p2 = p1 { }
  local p3 = p1 { }
  local p4 = p3 { }

  assertFalse( object:is_instanceOf(object) )

  assertTrue ( p0:is_instanceOf(object) )
  assertFalse( p0:is_instanceOf(p0) )
  assertFalse( p0:is_instanceOf(p1) )
  assertFalse( p0:is_instanceOf(p2) )
  assertFalse( p0:is_instanceOf(p3) )
  assertFalse( p0:is_instanceOf(p4) )

  assertTrue ( p1:is_instanceOf(object) )
  assertTrue ( p1:is_instanceOf(p0) )
  assertFalse( p1:is_instanceOf(p1) )
  assertFalse( p1:is_instanceOf(p2) )
  assertFalse( p1:is_instanceOf(p3) )
  assertFalse( p1:is_instanceOf(p4) )

  assertTrue ( p2:is_instanceOf(object) )
  assertTrue ( p2:is_instanceOf(p0) )
  assertTrue ( p2:is_instanceOf(p1) )
  assertFalse( p2:is_instanceOf(p2) )
  assertFalse( p2:is_instanceOf(p3) )
  assertFalse( p2:is_instanceOf(p4) )

  assertTrue ( p3:is_instanceOf(object) )
  assertTrue ( p3:is_instanceOf(p0) )
  assertTrue ( p3:is_instanceOf(p1) )
  assertFalse( p3:is_instanceOf(p2) )
  assertFalse( p3:is_instanceOf(p3) )
  assertFalse( p3:is_instanceOf(p4) )

  assertTrue ( p4:is_instanceOf(object) )
  assertTrue ( p4:is_instanceOf(p0) )
  assertTrue ( p4:is_instanceOf(p1) )
  assertFalse( p4:is_instanceOf(p2) )
  assertTrue ( p4:is_instanceOf(p3) )
  assertFalse( p4:is_instanceOf(p4) )

  assertFalse( is_instanceOf(0 , p0) )
  assertFalse( is_instanceOf('', p0) )
  assertFalse( is_instanceOf({}, p0) )
end

function TestLuaObjectErr:testSetFlag()
  local p0 = object 'p0' {}
  local msg = {
    "invalid argument #2 (number expected)"
  }

  for i=2,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set_flag, objectErr[i])
    assertErrorMsgContains( msg[1], p0.set_flag, p0, objectErr[i])
  end
end

function TestLuaObjectErr:testClearFlag()
  local p0 = object 'p0' {}
  local msg = {
    "invalid argument #2 (number expected)"
  }

  for i=2,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.clear_flag, objectErr[i])
    assertErrorMsgContains( msg[1], p0.clear_flag, p0, objectErr[i])
  end
end

function TestLuaObjectErr:testTestFlag()
  local p0 = object 'p0' {}
  local msg = {
    "invalid argument #2 (number expected)"
  }

  for i=2,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.test_flag, objectErr[i])
    assertErrorMsgContains( msg[1], p0.test_flag, p0, objectErr[i])
  end
end

function TestLuaObject:testSetGetFlags()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object {} :set_readonly()
  local p3 = object
  local p00 = p0 'p00' {}

  assertEquals(p0:set_flags(-1), p0)
  assertEquals(p1:set_flags(-1), p1)
  assertEquals(p2:set_flags(-1), p2)
  assertEquals(p3:set_flags(-1), p3)

  local tobit in MAD.operator
  assertEquals(p0:get_flags(), tobit(0xfffffffe)) -- object (class+        )
  assertEquals(p1:get_flags(), tobit(0xfffffffc)) -- object (     +        )
  assertEquals(p2:get_flags(), tobit(0xfffffffd)) -- object (     +readonly)
  assertEquals(p3:get_flags(), tobit(0xffffffff)) -- object (class+readonly)

  -- no inheritance
  p0:set_flags(-1)
  assertEquals(p00:get_flags(), 0)
  p00:set_flags(0)
  assertNotEquals(p0:get_flags(), 0)
  local p01 = p0 'p00' {}
  assertNotEquals(p01:get_flags(), p0:get_flags())
end

function TestLuaObjectErr:testSetGetFlags()
  local p0 = object 'p0' {}
  local msg = {
    "invalid argument #2 (number expected)"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set_flags, objectErr[i])
    assertErrorMsgContains(_msg[1], p0.get_flags, objectErr[i])
  end
  assertErrorMsgContains(msg[1], p0.set_flags, p0, nil)
  assertErrorMsgContains(msg[1], p0.set_flags, p0, true)
  assertErrorMsgContains(msg[1], p0.set_flags, p0, "")
  assertErrorMsgContains(msg[1], p0.set_flags, p0, {})
  assertErrorMsgContains(msg[1], p0.set_flags, p0, myFunc)
end

function TestLuaObject:testClearSetTestFlag()
  local p0 = object 'p0' {}
  local p1 = object {}
  local p2 = object
  local p00 = p0 'p00' {}

  assertFalse(p0:test_flag(2))
  assertFalse(p1:test_flag(2))
  assertFalse(p2:test_flag(2))
  for i=2, 31 do
    assertFalse(p0:test_flag(i))
    p0:set_flag(i)
    assertTrue(p0:test_flag(i))
    p0:clear_flag(i)
    assertFalse(p0:test_flag(i))
  end

  -- no inheritance
  p0:set_flag(2)
  assertFalse(p00:test_flag(2))
  p00:clear_flag(2)
  assertTrue(p0:test_flag(2))
  local p01 = p0 'p00' {}
  assertFalse(p01:test_flag(2))
end

function TestLuaObject:testValueSemantic()
  local p0 = object {}
  local p1 = p0 { x=3, y=2, z=function(s) return 2*s.y end}
  local p2 = p1 { x=2, y=function(s) return 3*s.x end }
  local p3 = p2 { x=function() return 5 end }
  local p4 = p3 { }
  local vs = {'x','y','z'}

  assertEquals   ( p0 , {} )
  assertEquals   ( p1:get(vs)     , { x=3, y=2, z=4 } )
  assertNotEquals( p1:get(vs,true), { x=3, y=2, z=4 } )
  assertNotEquals( p1:get(vs)     , { x=3, y=2 } )
  assertEquals   ( p2:get(vs)     , { x=2, y=6, z=12 } )
  assertNotEquals( p2:get(vs,true), { x=2, y=6, z=12 } )
  assertNotEquals( p2:get(vs)     , { x=2, y=6 } )
  assertNotEquals( p2:get(vs)     , { x=3, y=2 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=15, z=30 } )
  assertNotEquals( p3:get(vs,true), { x=5, y=15, z=30 } )
  assertNotEquals( p3:get(vs)     , { x=5, y=15 } )
  assertNotEquals( p3:get(vs)     , { x=2, y=6 } )
  assertEquals   ( p4:get(vs)     , { x=5, y=15, z=30 } )
  assertNotEquals( p4:get(vs,true), { x=5, y=15, z=30 } )
  assertNotEquals( p4:get(vs)     , { x=5, y=15 } )

  p1.z = 6
  assertEquals   ( p1:get(vs)     , { x=3, y=2, z=6 } )
  assertEquals   ( p1:get(vs,true), { x=3, y=2, z=6 } )
  assertNotEquals( p1:get(vs)     , { x=3, y=2 } )
  assertEquals   ( p2:get(vs)     , { x=2, y=6, z=6 } )
  assertNotEquals( p2:get(vs,true), { x=2, y=6, z=6 } )
  assertNotEquals( p2:get(vs)     , { x=2, y=6 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=15, z=6 } )
  assertNotEquals( p3:get(vs,true), { x=5, y=15, z=6 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=15, z=6 } )
  assertEquals   ( p4:get(vs)     , { x=5, y=15, z=6 } )
  assertNotEquals( p4:get(vs,true), { x=5, y=15, z=6 } )
  assertNotEquals( p4:get(vs)     , { x=5, y=15 } )

  p2.y = 5
  assertEquals   ( p2:get(vs)     , { x=2, y=5, z=6 } )
  assertEquals   ( p2:get(vs,true), { x=2, y=5, z=6 } )
  assertNotEquals( p2:get(vs)     , { x=2, y=5 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=5, z=6 } )
  assertNotEquals( p3:get(vs,true), { x=5, y=5, z=6 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=5, z=6 } )
  assertEquals   ( p4:get(vs)     , { x=5, y=5, z=6 } )
  assertNotEquals( p4:get(vs,true), { x=5, y=5, z=6 } )
  assertNotEquals( p4:get(vs)     , { x=5, y=5 } )

  p3.x = 3
  assertEquals   ( p3:get(vs)     , { x=3, y=5, z=6 } )
  assertEquals   ( p3:get(vs,true), { x=3, y=5, z=6 } )
  assertNotEquals( p3:get(vs)     , { x=3, y=5 } )
  assertEquals   ( p4:get(vs)     , { x=3, y=5, z=6 } )
  assertEquals   ( p4:get(vs,true), { x=3, y=5, z=6 } )
  assertNotEquals( p4:get(vs)     , { x=3, y=5 } )
end

function TestLuaObject:testArrayValueSemantic()
  local p0 = object {}
  local p1 = p0 { x=3, y=2, z=function(s) return { x=3*s.x, y=2*s.y } end }
  local p2 = p1 { x=2, y=function(s) return 2*s.x end }
  local p3 = p2 { x=function() return 5 end }
  local p4 = p3 {}
  local vs = {'x','y','z'}

  assertEquals   ( p0 , {} )
  assertEquals   ( p1:get(vs)     , { x=3, y=2, z={x=9, y=4} } )
  assertNotEquals( p1:get(vs,true), { x=3, y=2, z={x=9, y=4} } )
  assertNotEquals( p1:get(vs)     , { x=3, y=2 } )
  assertEquals   ( p2:get(vs)     , { x=2, y=4, z={x=6,y=8} } )
  assertNotEquals( p2:get(vs,true), { x=2, y=6, z={x=6,y=8} } )
  assertNotEquals( p2:get(vs)     , { x=2, y=6 } )
  assertNotEquals( p2:get(vs)     , { x=3, y=2 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=10, z={x=15,y=20} } )
  assertNotEquals( p3:get(vs,true), { x=5, y=10, z={x=15,y=20} } )
  assertNotEquals( p3:get(vs)     , { x=5, y=15 } )
  assertNotEquals( p3:get(vs)     , { x=2, y=6 } )
  assertEquals   ( p4:get(vs)     , { x=5, y=10, z={x=15,y=20} } )
  assertNotEquals( p4:get(vs,true), { x=5, y=15, z={x=15,y=20} } )
  assertNotEquals( p4:get(vs)     , { x=5, y=15 } )

  p1:set { x=function() return 7 end }
  assertEquals   ( p1:get(vs)     , { x=7, y=2, z={x=21,y=4} } )
  assertNotEquals( p1:get(vs,true), { x=7, y=2, z={x=21,y=4} } )
  assertNotEquals( p1:get(vs)     , { x=7, y=2 } )
  assertNotEquals( p1:get(vs)     , { y=2 } )
  assertEquals   ( p2:get(vs)     , { x=2, y=4, z={x=6,y=8} } )
  assertNotEquals( p2:get(vs,true), { x=2, y=4, z={x=6,y=8} } )
  assertNotEquals( p2:get(vs)     , { x=2, y=6 } )
  assertEquals   ( p3:get(vs)     , { x=5, y=10, z={x=15,y=20} } )
  assertNotEquals( p3:get(vs,true), { x=5, y=15, z=6 } )
  assertNotEquals( p3:get(vs)     , { x=5, y=15 } )
  assertEquals   ( p4:get(vs)     , { x=5, y=10, z={x=15,y=20} } )
  assertNotEquals( p4:get(vs,true), { x=5, y=15, z={x=15,y=20} } )
  assertNotEquals( p4:get(vs)     , { x=5, y=15 } )
end

function TestLuaObject:testSpecialVariable()
  local p0 = object 'p0' {}
  local p1 = p0 { x=3, y=function(s) return 2*s.x end, z=function(s) return { x=3*s.x, y=2*s.y } end }
  local p2 = p1 { x=2, y=function(s) return 4*s.x end }
  local p3 = p2 { x=function() return 5 end }
  local p4 = p3 {}

  assertTrue     ( p0.parent == object )
  assertTrue     ( p1.parent == p0 )
  assertTrue     ( p2.parent == p1 )
  assertTrue     ( p3.parent == p2 )
  assertTrue     ( p4.parent == p3 )
  assertTrue     ( p0.parent == p0.__par )
  assertTrue     ( p1.parent == p1.__par )
  assertTrue     ( p2.parent == p2.__par )
  assertTrue     ( p3.parent == p3.__par )
  assertTrue     ( p4.parent == p4.__par )

  assertEquals   ( p0.name , 'p0' )
  assertEquals   ( p1.name , 'p0' )
  assertEquals   ( p0      , { __id='p0' } )
  assertEquals   ( p1.__id , 'p0' )

  assertEquals   ( p0.parent.name , 'object' )
  assertTrue     ( p0.parent:is_readonly() )
  assertFalse    ( p0:is_readonly() )

  assertEquals   ( p1:var_raw'x'      , p1.x )
  assertEquals   ( p1:var_raw'y'(p1)  , p1.y )
  assertNotEquals( p1:var_raw'y'(p2)  , p1.y )
  assertEquals   ( p1:var_raw'z'(p1).x, p1.z.x )
  assertEquals   ( p1:var_raw'z'(p1).y, p1.z.y )
  assertNotEquals( p1:var_raw'z'(p2).y, p1.z.y )

  assertEquals   ( p2:var_raw'x'      , p2.x )
  assertEquals   ( p2:var_raw'y'(p2)  , p2.y )
  assertNotEquals( p2:var_raw'y'(p1)  , p2.y )
  assertEquals   ( p2:var_raw'z'(p2).x, p2.z.x )
  assertEquals   ( p2:var_raw'z'(p2).y, p2.z.y )
  assertNotEquals( p2:var_raw'z'(p3).y, p2.z.y )

  assertEquals   ( p3:var_raw'x'(p3)  , p3.x )
  assertEquals   ( p3:var_raw'y'(p3)  , p3.y )
  assertNotEquals( p3:var_raw'y'(p2)  , p3.y )
  assertEquals   ( p3:var_raw'z'(p3).x, p3.z.x )
  assertEquals   ( p3:var_raw'z'(p3).y, p3.z.y )
  assertEquals   ( p3:var_raw'z'(p4).y, p3.z.y )

  assertEquals   ( p4:var_raw'x'(p4)  , p4.x )
  assertEquals   ( p4:var_raw'y'(p4)  , p4.y )
  assertEquals   ( p4:var_raw'y'(p3)  , p4.y )
  assertEquals   ( p4:var_raw'z'(p4).x, p4.z.x )
  assertEquals   ( p4:var_raw'z'(p4).y, p4.z.y )

  assertEquals   ( p2.parent:var_raw'x'      , p1.x )
  assertEquals   ( p2.parent:var_raw'y'(p1)  , p1.y )
  assertEquals   ( p2.parent:var_raw'z'(p1).x, p1.z.x )
  assertEquals   ( p2.parent:var_raw'z'(p1).y, p1.z.y )
  assertNotEquals( p2.parent:var_raw'z'(p2).y, p1.z.y )

  assertEquals   ( p3.parent:var_raw'x'      , p2.x )
  assertEquals   ( p3.parent:var_raw'y'(p2)  , p2.y )
  assertEquals   ( p3.parent:var_raw'z'(p2).x, p2.z.x )
  assertEquals   ( p3.parent:var_raw'z'(p2).y, p2.z.y )
  assertNotEquals( p3.parent:var_raw'z'(p3).y, p2.z.y )

  assertEquals   ( p4.parent:var_raw'x'()    , p3.x )
  assertEquals   ( p4.parent:var_raw'y'(p3)  , p3.y )
  assertEquals   ( p4.parent:var_raw'z'(p3).x, p3.z.x )
  assertEquals   ( p4.parent:var_raw'z'(p3).y, p3.z.y )
  assertEquals   ( p4.parent:var_raw'z'(p4).y, p3.z.y )
end

function TestLuaObject:testIterators()
  local p0 = object 'p0' { 2, function() return 3 end, 4, x=1, y=2, z=function(s) return s.x*3 end }
  local p1 = p0 'p1' { 7, function() return 8 end, x=-1, y={} }
  local c

  assertEquals(#p0, 3)
  assertEquals(#p1, 2)

  -- bypass function evaluation, v may be a function but loops get same length
  c=0 for k,v in  pairs(p0) do c=c+1 assertEquals(p0:var_raw(k), v) end
  assertEquals(c , 7)
  c=0 for k,v in  pairs(p1) do c=c+1 assertEquals(p1:var_raw(k), v) end
  assertEquals(c , 5)
  c=0 for i,v in ipairs(p0) do c=c+1 assertEquals(p0:var_raw(i), v) end
  assertEquals(c , 3)
  c=0 for i,v in ipairs(p1) do c=c+1 assertEquals(p1:var_raw(i), v) end
  assertEquals(c , 2)
  assertEquals(p0[1], 2)
  assertTrue( is_function(p0[2]) )
  assertEquals(p0[3], 4)
  assertEquals(p1[1], 7)
  assertTrue( is_function(p1[2]) )
end

function TestLuaObject:testGetVariables()
  local f = function() return 4 end
  local p0 = object 'p0' { x=1, y=2, z=function() return 3 end, z2=f}
  local p1 = p0 'p1' { x=-1, y={} }
  local vs = {'name', 'x','y','z'}

  assertEquals ( p0:get(vs), { name='p0', x=1 , y=2 , z=3 } )
  assertEquals ( p1:get(vs), { name='p1', x=-1, y={}, z=3 } )

  assertEquals ( p0.get(p0,vs), { name='p0', x=1 , y=2 , z=3 } )
  assertEquals ( p1.get(p1,vs), { name='p1', x=-1, y={}, z=3 } )

  assertEquals( p0:get(vs    , false).z , 3 )
  assertEquals( p0:get(vs    , nil  ).z , 3 )
  assertEquals( p0:get({'z2'}       ).z2, 4 )
  assertEquals( p0:get({'z2'}, true ).z2, f )

  assertTrue  ( is_function(p0:get(vs, true ).name) )
  assertTrue  ( is_function(p0:get(vs, true ).z) )
  assertFalse ( is_function(p0:get(vs, nil  ).z) )

  assertTrue  ( is_function(p1:get(vs, true ).name) )
  assertTrue  ( is_function(p1:get(vs, true ).z) )
  assertTrue  ( is_table   (p1:get(vs, true ).y) )
  assertTrue  ( is_table   (p1:get(vs, true ).y) )
end

function TestLuaObjectErr:testGetVariables()
  local p0 = object 'p0' { }
  local msg = {
    "invalid argument #2 (iterable expected)"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.get, objectErr[i])
  end
  assertErrorMsgContains(msg[1], p0.get, p0, 0)
  assertErrorMsgContains(msg[1], p0.get, p0, true)
  assertErrorMsgContains(msg[1], p0.get, p0, "")
  assertErrorMsgContains(msg[1], p0.get, p0, myFunc)
end

function TestLuaObject:testSetVariables()
  local p0 = object 'p0' {}
  local p1 = p0 'p1' {}
  local vs = {'name', 'x','y','z'}

  p0:set { x=1, y=2, z=function() return 3 end }
  p1:set { x=-1, y={} }

  assertEquals ( p0:get(vs), { name='p0', x=1 , y=2 , z=3 } )
  assertEquals ( p1:get(vs), { name='p1', x=-1, y={}, z=3 } )

  assertEquals ( p0.get(p0,vs), { name='p0', x=1 , y=2 , z=3 } )
  assertEquals ( p1.get(p1,vs), { name='p1', x=-1, y={}, z=3 } )

  assertEquals( p0:get(vs, false).z, 3 )
  assertEquals( p0:get(vs, nil  ).z, 3 )

  assertTrue  ( is_function(p0:get(vs, true ).name) )
  assertTrue  ( is_function(p0:get(vs, true ).z) )
  assertFalse ( is_function(p0:get(vs, nil  ).z) )

  assertTrue  ( is_function(p1:get(vs, true ).name) )
  assertTrue  ( is_function(p1:get(vs, true ).z) )
  assertTrue  ( is_table   (p1:get(vs, true ).y) )
  assertTrue  ( is_table   (p1:get(vs, true ).y) )

  assertEquals( p1:set({x=-2}        ).x, -2 )
  assertEquals( p1:set({x=-3}, true  ).x, -3 )
  assertEquals( p1:set({x=-4}, nil   ).x, -4 )
  assertEquals( p1:set({x=-4}, nil   ).x, -4 )
  assertEquals( p1:set({z=-5}, false ).z, -5 )

  p1:set { z=5 }
  assertEquals( p1:get(vs, true).z, 5 )
end

function TestLuaObjectErr:testSetVariables()
  local p0 = object 'p0' { x=1 }
  local p1 = p0 'p1' { y={} }
  local p2 = p0 'p2' {}:set_readonly()
  local vs = {'name', 'x','y','z'}
  local msg = {
    "invalid argument #2 (mappable expected)",
    "cannot override variable"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.set, object, {})
  assertErrorMsgContains(_msg[3], p2.set, p2, {})

  assertErrorMsgContains(msg[1], p0.set, p0, 0)
  assertErrorMsgContains(msg[1], p0.set, p0, nil)
  assertErrorMsgContains(msg[1], p0.set, p0, true)
  assertErrorMsgContains(msg[1], p0.set, p0, "")
  assertErrorMsgContains(msg[1], p0.set, p0, myFunc)

  assertErrorMsgContains(msg[2], p0.set, p0, {x=2}, false)
  assertErrorMsgContains(msg[2], p1.set, p1, {y=2}, false)
end

function TestLuaObject:testWrapVariables()
  local f = function() return 3 end
  local p0 = object 'p0' { x=1, y=2, z=f }
  local p1 = p0 'p1' { x=-1, y={} }
  local vw = {'name', 'x','y','z'}

  -- very simple example
  p0:wrap_variables{x=function(e) return e()+2 end, y =function(e) return e()^3 end}
  assertEquals(p0:get{'x', 'y'}, {x=3, y=8})

  -- inheritance
  p1:wrap_variables{z=function(e) return function() return e end end}
  assertEquals(p1.z  , f)
  assertEquals(p1.z(), 3)

  -- keep functor semantic
  p0:set_methods{a=function(s) return s.name end}
  p0:wrap_variables{a=function(prev) return function(s) return prev(s.parent)end end}
  assertEquals(p0:a(), "object")

  -- concrete example: change angle
  local eps = 2.2204460492503130e-16
  local ksb = 0.85
  local rbend = object 'rbend' {}
  local r1 = rbend 'r1' {
    angle =function() return ksb end,
    h=3,
    length=function(s) return 2*s.h end
  }
  local abs = function(x) return is_number(x) and math.abs(x) or x:abs() end
  local sinc = function(x) return abs(x)<1e-10 and 1 or math.sin(x)/x end
  r1:wrap_variables {
    length = function(l) return function(s) return l(s)/sinc(s.angle)end end
  }
  assertAlmostEquals(r1.length, 6.78841077859289577, 8*eps)
end

function TestLuaObjectErr:testWrapVariables()
  local p0 = object 'p0' { y=1 }
  local p1 = p0 'p1' {}:set_readonly()
  local msg = {
    "invalid argument #2 (mappable expected)",
    "invalid variable (nil value)",
    "invalid wrapper (callable expected)",
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.wrap_variables, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.wrap_variables, object, {})
  assertErrorMsgContains(_msg[3], p1.wrap_variables, p1, {})

  assertErrorMsgContains(msg[1], p0.wrap_variables, p0, 0)
  assertErrorMsgContains(msg[1], p0.wrap_variables, p0, nil)
  assertErrorMsgContains(msg[1], p0.wrap_variables, p0, true)
  assertErrorMsgContains(msg[1], p0.wrap_variables, p0, "")
  assertErrorMsgContains(msg[1], p0.wrap_variables, p0, myFunc)

  assertErrorMsgContains(msg[2], p0.wrap_variables, p0, {z=0})
  assertErrorMsgContains(msg[2], p0.wrap_variables, p0, {z=true})
  assertErrorMsgContains(msg[2], p0.wrap_variables, p0, {z=""})
  assertErrorMsgContains(msg[2], p0.wrap_variables, p0, {z={}})
  assertErrorMsgContains(msg[2], p0.wrap_variables, p0, {z2=function() return 2 end})

  assertErrorMsgContains(msg[3], p0.wrap_variables, p0, {y=0})
  assertErrorMsgContains(msg[3], p0.wrap_variables, p0, {y=true})
  assertErrorMsgContains(msg[3], p0.wrap_variables, p0, {y=""})
  assertErrorMsgContains(msg[3], p0.wrap_variables, p0, {y={}})
end

function TestLuaObject:testSetFunction()
  local p0 = object 'p0' { z=function() return 3 end }
  local p1 = p0 'p1' {}

  p0:set_methods { x=function() return 2 end, y=function(s,n) return s.z*n end }

  assertFalse ( is_function(p0.x) )
  assertFalse ( is_function(p0.y) )
  assertFalse ( is_function(p0.z) )
  assertEquals( p0.z, 3 )
  assertFalse ( is_function(p0:get{'z'}.z) )
  assertTrue  ( is_function(p0:get({'z'},true).z) )
  assertTrue  ( is_function(p0:var_raw'z') )

  assertFalse ( is_function(p1.x) )
  assertFalse ( is_function(p1.y) )
  assertFalse ( is_function(p1.z) )
  assertEquals( p1.x(), 2)
  assertEquals( p1:y(3), 9)

  assertTrue  ( is_functor(p1.x) )
  assertTrue  ( is_functor(p1.y) )
  assertTrue  ( is_number(p1.z) )

  p1:set_methods({ x=function() return function() return 2 end end, y =function(s) return function(n) return s.z*n end end })

  assertTrue  ( is_functor (p1.x)   )
  assertTrue  ( is_function(p1.x()) )
  assertEquals( p1.x()(), 2)

  assertTrue  ( is_functor   (p1.y) )
  assertTrue  ( is_function(p1:y()) )
  assertEquals( p1:y()(3), 9)

  p1.y = function(s) return function(n) return s.z*n end end
  assertFalse ( is_functor   (p1.y) )
  assertTrue  ( is_function(p1.y) )
  assertFalse ( is_function(p1.y(3)) )
  assertEquals( p1.y(3), 9)
end

function TestLuaObjectErr:testSetMethods()
  local p0 = object 'p0' { x=function() return 2 end }
  local p1 = p0 'p1' { y=function() return 2 end }
  local p2 = p0 'p2' {}:set_readonly()
  local msg = {
    "invalid argument #2 (mappable expected)",
    "invalid key (function name expected)",
    "invalid value (callable expected)",
    "cannot override function"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set_methods, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.set_methods, object, {})
  assertErrorMsgContains(_msg[3], p2.set_methods, p2, {})

  assertErrorMsgContains(msg[1], p0.set_methods, p0, 0)
  assertErrorMsgContains(msg[1], p0.set_methods, p0, nil)
  assertErrorMsgContains(msg[1], p0.set_methods, p0, true)
  assertErrorMsgContains(msg[1], p0.set_methods, p0, "")
  assertErrorMsgContains(msg[1], p0.set_methods, p0, myFunc)

  assertErrorMsgContains(msg[2], p0.set_methods, p0, {myFunc, 2, {3}})

  assertErrorMsgContains(msg[3], p0.set_methods, p0, {x=0})
  assertErrorMsgContains(msg[3], p0.set_methods, p0, {x=true})
  assertErrorMsgContains(msg[3], p0.set_methods, p0, {x=""})
  assertErrorMsgContains(msg[3], p0.set_methods, p0, {x={}})

  assertErrorMsgContains(msg[4], p0.set_methods, p0, {x=function() return 3 end}, false)
  assertErrorMsgContains(msg[4], p1.set_methods, p1, {y=function() return 3 end}, false)
end

function TestLuaObject:testSetMetamethod()
  local p0 = object 'p0' { 1, 2, z=function() return 3 end }
  local p1 = p0 'p1' {}
  local tostr = function(s)
      local str = ''
      for k,v in pairs(s) do str = str .. tostring(k) .. ', ' end
      return str .. '#=' .. tostring(#s)
    end

  local p00 = object 'p00' { 1, 2, z=function() return 3 end }
  -- clone metatable shared with object
  p00:set_metamethods({ __tostring = tostr }, true)
  assertEquals      (tostring(p00), '1, 2, z, __id, #=2') -- tostring -> tostr
  assertNotEquals   (tostring(p1), '__id, #=2')           -- builtin tostring
  assertStrContains (tostring(p1), "object: 'p1'")        -- builtin tostring

  -- p1 child of p0 and not yet a class
  -- clone metatable shared with object
  p1:set_metamethods({ __tostring = tostr }, true)
  assertEquals      (tostring(p1), '__id, #=2')        -- tostring -> tostr

  local p01 = p00 'p01' {} -- fresh p1
  p01:set_metamethods({ __tostring = tostr }, true)    -- clone need override
  assertEquals      (tostring(p01), '__id, #=2')       -- tostring -> tostr
  p01:set_methods { x=function() return 2 end, y =function(s) return function(n) return s.z*n end end}
  p01.z =function() return 3 end
  assertEquals      (tostring(p01), 'y, x, __id, z, #=2') -- tostring -> tostr

  -- example of the help
  local obj1 = object 'obj1' { e = 3 }
  obj1:set_metamethods({ __pow = function(s,p) return s.e^p end })
  assertEquals(obj1^3, 27)
  obj1:set_metamethods({ __tostring =function(s) return s.name .. " has e = " .. s.e end}, true)
  assertEquals(tostring(obj1), "obj1 has e = 3")
end

function TestLuaObjectErr:testSetMetamethod()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }
  local p2 = p0 'p2' { }:set_readonly()
  local msg = {
    "invalid argument #2 (mappable expected)",
    "invalid metatable write access (unexpected class)",
    "invalid key (metamethod expected)",
    "cannot override metamethod"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.set_metamethods, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.set_metamethods, object, {})
  assertErrorMsgContains(_msg[3], p2.set_metamethods, p2, {})

  assertErrorMsgContains(msg[1], p0.set_metamethods, p0, 0)
  assertErrorMsgContains(msg[1], p0.set_metamethods, p0, nil)
  assertErrorMsgContains(msg[1], p0.set_metamethods, p0, true)
  assertErrorMsgContains(msg[1], p0.set_metamethods, p0, "")
  assertErrorMsgContains(msg[1], p0.set_metamethods, p0, myFunc)

  -- p1 created means p0 is a class and cannot modify its metatable
  assertErrorMsgContains(msg[2], p0.set_metamethods, p0, {})

  assertErrorMsgContains(msg[3], p0.set_metamethods, p1, {x=1})
  assertErrorMsgContains(msg[3], p1.set_metamethods, p1, {y=''})

  assertErrorMsgContains(msg[4], p1.set_metamethods, p1, {__index=false}, false)
end

function TestLuaObject:testGetVarKeys()
  local p0 = object 'p0' { x=1, y=2, z=function() return 3 end}
  local p1 = p0 'p1' { x=-1, y={} }
  local r  = {"parent", "first_free_flag", "name"}
  local r2 = {"first_free_flag", "name", "parent", 'x','y','z'}

  local t = object:get_varkeys() assertEquals ( t, r ) -- object not excluded
  t = p0:get_varkeys()             table.sort(t) assertEquals (t, r2)
  t = p0:get_varkeys(object)       table.sort(t) assertEquals (t, {'x','y','z'})
  t = p0:get_varkeys(object.__par) table.sort(t) assertEquals (t, r2)

  t = p1:get_varkeys()             table.sort(t) assertEquals (t, r2)
  t = p1:get_varkeys(p0.__par)     table.sort(t) assertEquals (t, {'x','y','z'})
  t = p1:get_varkeys(p0)           table.sort(t) assertEquals (t, {'x','y'})
  t = p1:get_varkeys(p1.__par)     table.sort(t) assertEquals (t, {'x','y'})

end

function TestLuaObjectErr:testGetVarKeys()
  local p0 = object 'p0' { x=1, y=2, z=function() return 3 end }
  local p1 = p0 'p1' { }
  local msg = {
    "invalid argument #2 (parent of argument #1 expected)"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.get_varkeys, objectErr[i])
  end
  assertErrorMsgContains(_msg[2], p0.get_varkeys, p0, 0)
  assertErrorMsgContains(_msg[2], p0.get_varkeys, p0, true)
  assertErrorMsgContains(_msg[2], p0.get_varkeys, p0, "")
  assertErrorMsgContains(_msg[2], p0.get_varkeys, p0, {})
  assertErrorMsgContains(_msg[2], p0.get_varkeys, p0, myFunc)

  assertErrorMsgContains(msg[1], p0.get_varkeys, p0, p1)
end

function TestLuaObject:testInsert()
  local p0 = object 'p0' { 1, 2, 3 }
  local p1 = p0 'p1' { }

  assertEquals(p0:insert(4,4), p0)
  assertEquals(p0[4], 4)
  assertEquals(p1[4], 4)

  assertEquals(p1:insert(5,5), p1)
  assertEquals(p1[5], 5)

  p0:insert(2,"test")
  assertEquals(p0[1], 1)
  assertEquals(p0[2], "test")
  assertEquals(p0[3], 2)
  assertEquals(p0[4], 3)

  p0:insert(1000,"test")
  assertEquals(p0[1000], "test")
end

function TestLuaObjectErr:testInsert()
  local p0 = object 'p0' { 1, 2, 3 }
  local p1 = p0 'p1' {}:set_readonly()
  local msg = {
    "bad argument #2 to 'insert' (number expected, got nil)",
    "bad argument #2 to 'insert' (number expected, got boolean)",
    "bad argument #2 to 'insert' (number expected, got string)",
    "bad argument #2 to 'insert' (number expected, got table)",
    "bad argument #2 to 'insert' (number expected, got function)",
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.insert, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.insert, object, 1)
  assertErrorMsgContains(_msg[3], p1.insert, p1, 1)

  assertErrorMsgContains(msg[1], p0.insert, p0, nil)
  assertErrorMsgContains(msg[2], p0.insert, p0, true)
  assertErrorMsgContains(msg[3], p0.insert, p0, "")
  assertErrorMsgContains(msg[4], p0.insert, p0, {})
  assertErrorMsgContains(msg[5], p0.insert, p0, myFunc)
end

function TestLuaObject:testRemove()
  local p0 = object 'p0' { 1, 2, 3, 4 }
  local p1 = p0 'p1' { }

  assertEquals(p0:remove(2), 2)
  assertEquals(p0[1], 1)
  assertEquals(p0[2], 3)
  assertEquals(p0[3], 4)
  assertEquals(p1[1], 1)
  assertEquals(p1[2], 3)
  assertEquals(p1[3], 4)

  assertEquals(p1:remove(2), nil)
  assertEquals(p0[1], 1)
  assertEquals(p0[2], 3)
  assertEquals(p0[3], 4)
  assertEquals(p1[1], 1)
  assertEquals(p1[2], 3)
  assertEquals(p1[3], 4)
end

function TestLuaObjectErr:testRemove()
  local p0 = object 'p0' { 1, 2, 3, 4 }
  local p1 = p0 'p1' {}:set_readonly()
  local msg = {
    "bad argument #2 to '?' (number expected, got boolean)",
    "bad argument #2 to '?' (number expected, got string)",
    "bad argument #2 to '?' (number expected, got table)",
    "bad argument #2 to '?' (number expected, got function)",
  }

  assertErrorMsgContains(_msg[3], object.remove, object, 1)
  assertErrorMsgContains(_msg[3], p1.remove, p1, 1)

  assertErrorMsgContains(msg[1], p0.remove, p0, true)
  assertErrorMsgContains(msg[2], p0.remove, p0, "")
  assertErrorMsgContains(msg[3], p0.remove, p0, {})
  assertErrorMsgContains(msg[4], p0.remove, p0, myFunc)
end

function TestLuaObject:testBsearch()
  local p0 = object 'p0' { 1, 2, 2, 2, 3, 4, 5 }


  local res, resEq = {}, {}
  for i=1,7 do
    res[i]   = p0:bsearch(i, lt)
    resEq[i] = p0:bsearch(i, le)
  end

  assertEquals(res  , {1, 2, 5, 6, 7, 8, 8})
  assertEquals(resEq, {2, 5, 6, 7, 8, 8, 8})

  res, resEq = {}, {}
  p0:sort(gt)
  for i=1,7 do
    res[i]   = p0:bsearch(i, gt)
    resEq[i] = p0:bsearch(i, ge)
  end

  assertEquals(res  , {7, 4, 3, 2, 1, 1, 1})
  assertEquals(resEq, {8, 7, 4, 3, 2, 1, 1})
end

function TestLuaObjectErr:testBsearch()
  local p0 = object 'p0' { }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.bsearch, objectErr[i])
  end
end

function TestLuaObject:testSort()
  local p0 = object 'p0' { 1, 5, 2, -10, 3, 4, 2, 200 }
  local p1 = p0 'p1' { }

  p0:sort(lt)
  assertEquals(p0:get{1,2,3,4,5,6,7,8}, {-10,1,2,2,3,4,5,200})
  p0:sort(gt)
  assertEquals(p0:get{1,2,3,4,5,6,7,8}, {200,5,4,3,2,2,1,-10})
end

function TestLuaObjectErr:testSort()
  local p0 = object 'p0' { 1, 5, 2, -10, 3, 4, 2, 200 }
  local p1 = p0 'p1' {}:set_readonly()
  local msg = {
    "bad argument #2 to 'sort' (function expected, got boolean)",
    "bad argument #2 to 'sort' (function expected, got number)",
    "bad argument #2 to 'sort' (function expected, got table)",
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.sort, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.sort, object)
  assertErrorMsgContains(_msg[3], p1.sort, p1)

  assertErrorMsgContains(msg[1], p0.sort, p0, true)
  assertErrorMsgContains(msg[2], p0.sort, p0, 1)
  assertErrorMsgContains(msg[3], p0.sort, p0, {})
end

function TestLuaObject:testClearArray()
  local p0 = object 'p0' {
      1,   true,   "",    {},    myFunc,
    x=1, y=true, z="", z2={}, z3=myFunc
  }
  local p1 = p0 'p1' { }
  local name = {"x", "y", "z", "z2", "z3"}
  local tbl = {x=1, y=true, z="", z2={}, z3 =myFunc}

  -- no inheritance
  assertEquals(p1:clear_array(), p1)
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p1, 5)
  assertEquals(#p0, 5)

  assertEquals(p0:clear_array(), p0)
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p1, 0)
  assertEquals(#p0, 0)

  p0:set({1, 2, 3, 4, 5})
  p1:set({1, 2, 3})
  assertEquals(#p0, 5)
  assertEquals(#p1, 3)
  p1:clear_array()
  assertEquals(#p0, 5)
  assertEquals(#p1, 5)
end

function TestLuaObjectErr:testClearArray()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' {}:set_readonly()

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.clear_array, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.clear_array, object)
  assertErrorMsgContains(_msg[3], p1.clear_array, p1)
end

function TestLuaObject:testClearVariables()
  local p0 = object 'p0' {
      1,   true,   "",    {},    myFunc,
    x=1, y=true, z="", z2={}, z3=myFunc
  }
  local p1 = p0 'p1' { }
  local name = {"x", "y", "z", "z2", "z3"}
  local tbl = {x=1, y=true, z="", z2={}, z3 =myFunc}

  -- no inheritance
  assertEquals(p1:clear_variables(), p1)
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p1, 5)
  assertEquals(#p0, 5)

  assertEquals(p0:clear_variables(), p0)
  assertEquals(p0:get(name, true), {})
  assertEquals(p1:get(name, true), {})
  assertEquals(#p1, 5)
  assertEquals(#p0, 5)

  p0:set(tbl)
  p1:set({x=2, y=false})
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get({"x", "y"}, true), {x=2, y=false})
  assertEquals(#p0, 5)
  assertEquals(#p1, 5)
  p1:clear_variables()
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p0, 5)
  assertEquals(#p1, 5)
end

function TestLuaObjectErr:testClearVariables()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' {}:set_readonly()

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.clear_variables, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.clear_variables, object)
  assertErrorMsgContains(_msg[3], p1.clear_variables, p1)
end

function TestLuaObject:testClearAll()
  local p0 = object 'p0' {
      1,   true,   "",    {},    myFunc,
    x=1, y=true, z="", z2={}, z3=myFunc
  }
  local p1 = p0 'p1' { }
  local name = {"x", "y", "z", "z2", "z3"}
  local tbl = {x=1, y=true, z="", z2={}, z3 =myFunc}

  -- no inheritance
  assertEquals(p1:clear_all(), p1)
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p1, 5)
  assertEquals(#p0, 5)

  assertEquals(p0:clear_all(), p0)
  assertEquals(p0:get(name, true), {})
  assertEquals(p1:get(name, true), {})
  assertEquals(#p1, 0)
  assertEquals(#p0, 0)

  p0:set(tbl)
  p0:set({1, 2, 3, 4, 5})
  p1:set({x=2, y=false})
  p1:set({1, 2, 3})
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get({"x", "y"}, true), {x=2, y=false})
  assertEquals(#p0, 5)
  assertEquals(#p1, 3)
  p1:clear_all()
  assertEquals(p0:get(name, true), tbl)
  assertEquals(p1:get(name, true), tbl)
  assertEquals(#p0, 5)
  assertEquals(#p1, 5)
end

function TestLuaObjectErr:testClearAll()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }:set_readonly()

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.clear_all, objectErr[i])
  end
  assertErrorMsgContains(_msg[3], object.clear_all, object)
  assertErrorMsgContains(_msg[3], p1.clear_all, p1)
end

function TestLuaObject:testRawLen()
  local p0 = object 'p0' { 1, 2, 3, 4 }
  local p1 = p0 'p1' { x=true, y=true }
  local p2 = p0 'p2' { 1, 2, 3 }

  assertEquals(#p0, 4)
  assertEquals(p0:raw_len(), 4)

  assertEquals(#p1, 4)
  assertEquals(p1:raw_len(), 0)

  assertEquals(#p2, 3)
  assertEquals(p2:raw_len(), 3)
end

function TestLuaObject:testRawGetSet()
  local p0 = object 'p0' { 1, 2, 3, 4 }
  local p1 = p0 'p1' {42, x=true, y=function() return true end}
  local p2 = p0 'p2' {}:set_readonly()

  p2:raw_set("x", true)
  assertTrue(p2.x)

  assertEquals(p1:raw_get(1), 42)
  assertEquals(p1:raw_get("x"), true)
  assertTrue  (is_function(p1:raw_get("y")))
  assertFalse (is_function(p1.y))
  assertNil   (p1:raw_get(2))
  assertNil   (p1:raw_get(3))
end

function TestLuaObjectErr:testRawGetSet()
  local p0 = object 'p0' {}:set_readonly()
  local msg = {
    "forbidden write access to 'p0' (readonly object or variable)",
  }

  assertErrorMsgContains(msg[1], function(p0) p0.x=false end, p0)
end

function TestLuaObject:testVarRawVal()
  local p0 = object 'p0' { }

  -- eval
  assertFalse (is_function(p0:var_val("k", function() return 2 end)))
  assertEquals(p0:var_val("k", function() return 2 end), 2)

  -- eval with self
  p0.x = 3
  assertFalse (is_function(p0:var_val("k", function(s) return s.x*2 end)))
  assertEquals(p0:var_val("k", function(s) return s.x*2 end), 6)

  -- no eval key not a strings
  assertTrue  (is_function(p0:var_val(1     , function() return 2 end)))
  assertTrue  (is_function(p0:var_val(true  , function() return 2 end)))
  assertTrue  (is_function(p0:var_val({}    , function() return 2 end)))
  assertTrue  (is_function(p0:var_val(myFunc, function() return 2 end)))
  assertTrue  (is_function(p0:var_val(p0    , function() return 2 end)))

  -- no eval value not a function
  assertEquals(p0:var_val("k", 1   ), 1   )
  assertEquals(p0:var_val("k", ""  ), ""  )
  assertEquals(p0:var_val("k", true), true)
  assertEquals(p0:var_val("k", {}  ), {})
end

function TestLuaObject:testVarRawGet()
  local p0 = object 'p0' { 1, 2, 3, 4, z=function() return 5 end }
  local p1 = p0 'p1' {42, x=true, y=function() return true end }

  assertEquals(p1:var_raw(1), 42)
  assertEquals(p1:var_raw(2), 2)
  assertEquals(p1:var_raw(3), 3)
  assertEquals(p1:var_raw("x"), true)
  assertTrue  (is_function(p1:var_raw("y")))
  assertFalse (is_function(p1.y))
  assertTrue  (is_function(p1:var_raw("z")))
  assertEquals(p1.z, 5)

  assertEquals(p1:var_get(1), 42)
  assertEquals(p1:var_get(2), 2)
  assertEquals(p1:var_get(3), 3)
  assertEquals(p1:var_get("x"), true)
  assertFalse (is_function(p1: var_get("y")))
  assertTrue  (p1:var_get("y"))
  assertFalse (is_function(p1:var_get("z")))
  assertEquals(p1:var_get("z"), 5)
end

function TestLuaObjectErr:testEnv()
  local p0 = object 'p0' { }
  local msg = {
    "invalid argument #2 (not a function or < 1)",
    "invalid environment (already open)",
    "invalid environment (not open)"
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.open_env, objectErr[i])
    assertErrorMsgContains(_msg[1], p0.is_open_env, objectErr[i])
    assertErrorMsgContains(_msg[1], p0.reset_env, objectErr[i])
    assertErrorMsgContains(_msg[1], p0.close_env, objectErr[i])
  end
  assertErrorMsgContains(msg[1], p0.open_env, p0, 0)
  assertErrorMsgContains(msg[1], p0.open_env, p0, -1)
  assertErrorMsgContains(msg[1], p0.open_env, p0, true)
  assertErrorMsgContains(msg[1], p0.open_env, p0, "")
  assertErrorMsgContains(msg[1], p0.open_env, p0, {})

  p0:open_env()
  assertErrorMsgContains(msg[2], p0.open_env, p0)
  p0:close_env()

  assertErrorMsgContains(msg[3], p0.close_env, p0)
end

function TestLuaObject:testEnvSimple()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }

  local oCtx, nCtx = getfenv(1), {}
  local setfenv = setfenv

  local res = {a=2, b="test", c={a=2, b="test"}, d=3}

  setfenv(1, nCtx)
    p0:open_env()
      a = 2
      b = "test"
      c = {a=a, b=b}
      d = function(s) return s.a+1 end
    p0:close_env()
    e = 4
  setfenv(1, oCtx)

  assertEquals(nCtx, {e=4})
  assertEquals(p0:get{"a", "b", "c", "d", "e"}, res)
end

function TestLuaObject:testEnvInheritance()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }

  assertFalse(p0:is_open_env())
  assertFalse(p1:is_open_env())

  p0:open_env()
    assertTrue (p0:is_open_env())
    assertFalse(p1:is_open_env())
  p0:close_env()

  assertFalse(p0:is_open_env())
  assertFalse(p1:is_open_env())
end

function TestLuaObject:testEnvNested()
  local p0 = object 'p0' { }
  local p1 = p0 'p1' { }

  local oCtx, nCtx = getfenv(1), {}
  local setfenv = setfenv

  local names = {"a", "a0", "a1", "b0", "b"}

  setfenv(1, nCtx)
    p0:clear_all()
    p1:clear_all()
    a = true
    p0:open_env()
      a0 = true
      p1:open_env()
        a1= true
      p1:close_env()
      b0 = true
    p0:close_env()
    b = true
  setfenv(1, oCtx)

  assertEquals(nCtx, {a=true, b=true})
  assertEquals(p0:get(names), {a0=true, b0=true})
  assertEquals(p1:get(names), {a0=true, a1=true, b0=true})
  assertNil(p1:raw_get("a0"))
  assertTrue(p1:raw_get("a1"))
  assertNil(p1:raw_get("b0"))
end

function TestLuaObject:testEnvFunc()
  local p0 = object 'p0' { }
  local oCtx, nCtx = getfenv(1), {}
  local setfenv = setfenv

  setfenv(1, nCtx)
    local function func()
      a=1
      b=true
      assertTrue(p0:is_open_env())
    end

    p0:open_env(func)
    func()
  setfenv(1, oCtx)

  assertEquals(p0:get{"a","b"}, {a=1, b=true})
  assertEquals(nCtx, {})
end

function TestLuaObject:testEnvMultLvl()
  local p0 = object 'p0' { }
  local oCtx, nCtx = getfenv(1), {}
  local setfenv = setfenv

  setfenv(1, nCtx)
    local function func1()
      a = true

      local function func2()
        b = true

        local function func3()
          c = true
          p0:open_env(func1)
          assertTrue(p0:is_open_env())
          d = true
          assertEquals({c,d}, {true,true})
        end

        func3()
        e = true
        assertTrue(p0:is_open_env())
        assertEquals({b,e}, {true,true})
      end

      func2()
      f = true
      assertTrue(p0:is_open_env())
      assertEquals({a,f}, {nil,true})
      p0:close_env()
    end
    func1()
  setfenv(1, oCtx)

  assertEquals(p0:get{"a","b","c","d","e","f"}, {f=true})
end

function TestLuaObject:testEnvSelfRef()
  local p0 = object 'p0' { }

  p0:open_env()
    assertTrue  (p0:is_open_env())
    assertTrue  (is_object(p0.p0))
    assertEquals(p0.p0, p0)
  p0:close_env()

  assertFalse (p0:is_open_env())
  assertFalse (is_object(p0.p0))
  assertEquals(p0.p0)
end

function TestLuaObject:testEnvReset()
  local p0 = object 'p0' { }
  local error = error
  local func = function()
    p0:open_env()
    error("throw error")
  end

  local status, msg = pcall(func)
  assertEquals(status, false)
  assertNotNil(string.find(msg, "throw error"))
  p0:reset_env()
  assertFalse(p0:is_open_env())
  assertNil(p0.p0)
end

function TestLuaObjectErr:testDumpObj()
  local p0 = object 'p0' { x=1, y=2, z:=3 }
  local p1 = p0 'p1' {}
  local msg = {
    "invalid argument #2 (parent of argument #1 expected)",
    "invalid argument #3 (object expected)",
    "invalid argument #4 (string expected)",
  }

  for i=1,#objectErr+1 do
    assertErrorMsgContains(_msg[1], p0.dumpobj, objectErr[i])
  end
  assertErrorMsgContains(msg[1], p0.dumpobj, p0, '-', p1)

  assertErrorMsgContains(msg[2], p0.dumpobj, p0, '-', 0)
  assertErrorMsgContains(msg[2], p0.dumpobj, p0, '-', true)
  assertErrorMsgContains(msg[2], p0.dumpobj, p0, '-', {})
  assertErrorMsgContains(msg[2], p0.dumpobj, p0, '-', myFunc)

  assertErrorMsgContains(msg[3], p0.dumpobj, p0, '-', object, 0)
  assertErrorMsgContains(msg[3], p0.dumpobj, p0, '-', object, true)
  assertErrorMsgContains(msg[3], p0.dumpobj, p0, '-', object, {})
  assertErrorMsgContains(msg[3], p0.dumpobj, p0, '-', object, myFunc)
  assertErrorMsgContains(msg[3], p0.dumpobj, p0, '-', object, object)
end

function TestLuaObject:testDumpObj()
  local p0 = object 'p0' { x=1, y=2, z:=3 }
  local p1 = p0 'p1' {}
  local p2 = p0 'p2' { x=-1, y={}, z2="", z3=\s s.x, z4=p1}
  local str_p0 = [[
+ object: 'p0'
   y :  2
   x :  1
   z := 3]]
  local str_p1 = [[
+ object: 'p2'
   x :  -1
   y :  {}
   z4 :  object: 'p1'
   z2 : ''
   z3 := -1
   + object: 'p0'
      y :  2 (*)
      x :  1 (*)
      z := 3]]
  local str_p1_np0 = [[
+ object: 'p2'
   x :  -1
   y :  {}
   z4 :  object: 'p1'
   z2 : ''
   z3 := -1]]
  local str_p1_pattern = [[
+ object: 'p2'
   z4 :  object: 'p1'
   z2 : ''
   z3 := -1
   + object: 'p0'
      z := 3]]

  assertEquals(string.gsub(p0:dumpobj('-'    ), '%s0x%x+' , ''), str_p0)
  assertEquals(string.gsub(p2:dumpobj('-'    ), '%s0x%x+' , ''), str_p1)
  assertEquals(string.gsub(p2:dumpobj('-', p0), '%s0x%x+' , ''), str_p1_np0)
  assertEquals(string.gsub(p2:dumpobj('-', nil, "z[0-9]?"), '%s0x%x+', ''),
               str_p1_pattern)
end

-- examples test suite --------------------------------------------------------o
function TestLuaObject:testMetamethodForwarding()
  local msg = {
    "invalid argument #1 (forbidden access to 'ro_obj')",
  }

  local ro_obj = object {}
  local parent = ro_obj.parent

  ro_obj:set_methods {
    set_readonly = function(s,f)
      assert(s ~= ro_obj, msg[1])
      return parent.set_readonly(s,f)
    end
  }
  ro_obj:set_metamethods { __init = function(s) return s:set_readonly(true) end}
  assertErrorMsgContains(msg[1], ro_obj.set_readonly, ro_obj, true)
  assertFalse( ro_obj:is_readonly() )
  parent.set_readonly(ro_obj, true)
  assertTrue ( ro_obj:is_readonly() )

  local ro_chld = ro_obj {}
  assertTrue ( ro_chld:is_readonly() )
  assertTrue ( ro_chld:set_readonly(1):is_readonly() )
  assertFalse( ro_chld:set_readonly(false):is_readonly() )
end

function TestLuaObject:testMetamethodNotification()
  local p1 = object 'p1' { x=1, y=2  }
  local p2 = p1 'p2' { x=2, y=-1, z=0 }

  local function trace (fp, self, k, v)
--[[fp:write("object: '", self.name,
             "' is updated for key: '", tostring(k),
             "' with value: ")
    if type(v) == "string"
      then fp:write(": '", tostring(v), "'\n")
      else fp:write(":  ", tostring(v),  "\n") end
]]  end

  local function set_notification (self, file)
    local fp = file or io.stdout
    local mt = getmetatable(self)
    local nwidx = mt and rawget(mt, '__newindex')
    local mm = function (self, k, v)
      trace(fp, self, k, v) -- logging
      nwidx(    self, k, v) -- forward
    end
    self:set_metamethods({__newindex=mm}, true) -- override!
  end

  set_notification(p2) -- new metamethod created, metatable is cloned
  p2.x = 3 -- new behavior, notify about update

  local p3 = p2 'p3' { x=3  } -- new, inherit metatable
  p3.x = 4 -- new behavior, notify about update

  local p4 = p2 'p4' { x=4 } -- new, inherit metatable
  p4.x = 5 -- new behavior, notify about update
end

function TestLuaObject:testMetamethodCounting()
  local count = 0
  local set_counter = function(s) s:set_metamethods {
    __init = function(s) count=count+1 ; return s end
  } end

  local o0 = object 'o0' {}          set_counter(o0)
  local o1 = o0 'o1' { a = 2 }       assertEquals( count, 1 )
  local o2 = o1 'o2' { a = 2 }       assertEquals( count, 2 )
  local a = object 'a' { x = o2.a }  assertEquals( count, 2 )
end

-- performance test suite -----------------------------------------------------o

Test_LuaObject = {}

function Test_LuaObject:testPrimes()
  local Primes = object {}

  Primes:set_methods {
    isPrimeDivisible = function(s,c)
      for i=3, s.prime_count do
        if s.primes[i] * s.primes[i] > c then break end
        if c % s.primes[i] == 0 then return true end
      end
      return false
    end,

    addPrime = function(s,c)
      s.prime_count = s.prime_count + 1
      s.primes[s.prime_count] = c
    end,

    getPrimes = function(s,n)
      s.prime_count, s.primes = 3, { 1,2,3 }
      local c = 5
      while s.prime_count < n do
        if not s:isPrimeDivisible(c) then
          s:addPrime(c)
        end
        c = c + 2
      end
    end
  }

  local p = Primes {}
  local t0 = os.clock()
  p:getPrimes(2e5)
  local dt = os.clock() - t0
  assertEquals( p.primes[p.prime_count], 2750131 )
  assertAlmostEquals( dt , 0.5, 1 )
end

function Test_LuaObject:testDuplicates()

  local DupFinder = object {}

  DupFinder:set_methods {
    find_duplicates = function(s,res)
      for _,v in ipairs(s) do
        res[v] = res[v] and res[v]+1 or 1
      end
      for _,v in ipairs(s) do
        if res[v] and res[v] > 1 then
          res[#res+1] = v
        end
        res[v] = nil
      end
    end,

    clear = function(s)
      for i=1,#s do s[i]=nil end
      return s
    end
  }

  local inp = DupFinder {'b','a','c','c','e','a','c','d','c','d'}
  local out = DupFinder {'a','c','d'}
  local res = DupFinder {}

  local t0 = os.clock()
  for i=1,5e5 do inp:find_duplicates(res:clear()) end
  local dt = os.clock() - t0
  assertEquals( res, out )
  assertAlmostEquals( dt , 0.5, 1 ) -- fails very often
end

function Test_LuaObject:testDuplicates2()
  local DupFinder = object {}
  local _len = {}

  DupFinder:set_methods {
    find_duplicates = function(s,res)
      for _,v in ipairs(s) do
        res[v] = (res[v] or 0) + 1
      end
      for _,v in ipairs(s) do
        if res[v] and res[v] > 1 then
          local len = res[_len]+1
          res[len], res[_len] = v, len
        end
        res[v] = nil
      end
    end,

    clear = function(s)
      for i=1,s[_len] do s[i]=nil end
      s[_len] = 0
      return s
    end
  }

  local inp = DupFinder {'b','a','c','c','e','a','c','d','c','d'}
  local out = DupFinder {'a','c','d'}
  local res = DupFinder { [_len]=0 }

  local t0 = os.clock()
  for i=1,5e5 do inp:find_duplicates(res:clear()) end
  local dt = os.clock() - t0
  assertEquals( res, out )
  assertAlmostEquals( dt , 0.5, 1 )
end

function Test_LuaObject:testLinkedList()
  local List = object {}
  local nxt = {}

  local function generate(n)
    local t = List {x=1}
    for j=1,n do t = List {[nxt]=t} end
    return t
  end

  local function find(t,k)
    if t[k] ~= nil then return t[k] end
    return find(t[nxt],k)
  end

  local l, s, n = generate(10), 0, 1e6
  local t0 = os.clock()
  for i=1,n do s = s + find(l, 'x') end
  local dt = os.clock() - t0
  assertEquals( s, n )
  assertAlmostEquals( dt, 0.5, 1 )
end

function Test_LuaObject:testSlides()
  local Point = object 'Point' { }
  local p1 = Point 'p1' { x=3, y=2, z=1 }
  local p2 = p1 'p2' { x=2, y=1 }
  local p3 = p2 'p3' { x=1 }
  local p4 = p3 'p4' { }
  local s, t0, dt

  s = 0
  t0 = os.clock()
  for i=1,5e8 do
    s = p1.z + p2.z + p3.z + p4.z + s
  end
  dt = os.clock() - t0
  assertEquals( s, 2000000000 )
  assertAlmostEquals( dt, 0.5, 1 )

  p1.getz = function(s) return s.z end

  s=0
  t0 = os.clock()
  for i=1,5e8 do
    s = p1.getz + p2.getz + p3.getz + p4.getz + s
  end
  dt = os.clock() - t0
  assertEquals( s, 2000000000 )
  assertAlmostEquals( dt, 0.5, 1 )

  p1.getz = nil
  p1:set_methods { getz = function(s) return s.z end }

  s=0
  t0 = os.clock()
  for i=1,5e8 do
    s = p1:getz() + p2:getz() + p3:getz() + p4:getz() + s
  end
  dt = os.clock() - t0
  assertEquals( s, 2000000000 )
  assertAlmostEquals( dt, 0.5, 1 )
end

-- end ------------------------------------------------------------------------o

-- run as a standalone test with luajit
if MAD == nil then
  os.exit( utest.LuaUnit.run() )
end
