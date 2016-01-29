--[=[
 o----------------------------------------------------------------------------o
 |
 | TFS table module
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
  tfstable -- TFS table

SYNOPSIS
  table = require "tfstable"
  my_tab = table 'mytab' { column_name_list... }

DESCRIPTION
  The module table creates TFS tables used by MAD. The columns can be
  accessed by name or by index and their content will grow automatically.
  If a name of a column is enclosed into a list, then its elements can be
  used to access the row.

RETURN VALUE
  The TFS table.

EXAMPLES
  table = require "tfstable"
  tab = table 'survey' { {'name'}, 'x', 'y', 'z', 'phi', 'theta', 'rho' }
  tab:add{ 'drift', 0.1, 0.2, 0.5, 0, 0, 0 }
  tab:add{ name='mq', x=0.2, y=0.4, z=1, phi=0, theta=0, rho=0 }
  tab:write()         -- equivalent to tab:write"survey.tfs"
  print(tab.x[2])     -- x of 'mq'
  print(tab.mq.x)

SEE ALSO
  sequence, element, beam
]]

-- requires ------------------------------------------------------------------o

local object = require 'object'
local utils  = require 'utils'

-- globals -------------------------------------------------------------------o

local _G_mt = getmetatable(_G)
setmetatable(_G, {
  __index = function(_, n)
    error("attempt to read to undeclared variable " .. n, 2)
  end,
  __newindex = function(_, n)
    error("attempt to write to undeclared variable " .. n, 2)
  end
})

-- locals --------------------------------------------------------------------o

local type, setmetatable = type, setmetatable
local pairs, ipairs = pairs, ipairs
local is_list = utils.is_list

-- metatable for the root of all tables
local MT = object {} 

 -- make the module the root of all tables
MT (M)
M.name = 'table'
M.kind = 'table'
M.is_table = true

-- functions -----------------------------------------------------------------o

-- special functions for [row_name][col_name] 'lazy' access and error handling
-- normal access is [col_idx|col_name][row_idx]

local function invalid_row(self, key)
  error("invalid table access for key '" .. key .. "', invalid row name")
end

local function invalid_col(self, key)
  error("invalid table access for key '" .. key .. "', invalid column name")
end

local function invalid_idx(self, key)
  error("invalid table access for key '" .. key .. "', invalid column index")
end

local function invalid_cnt(self, key)
  error("invalid table access for key '" .. key .. "', missing or invalid count number")
end

local function invalid_val(self, key)
  error("invalid table access for key '" .. key .. "', value is a table")
end

local function invalid_tbl(self, key)
  error("invalid table access for key '" .. key .. "', value is not a table")
end

local function invalid_len(self, key)
  error("invalid table access for row '" .. key .. "', invalid number of values")
end

local function invalid_key(self, key)
  error("invalid table access for row '" .. key .. "', invalid column name")
end

-- set values of existing row

local function add_refkey(self, idx)
  local key = self._refcol[idx]
  local col = self._refkey
  local ref = col[key]

      if      ref  ==  nil     then col[key]    = idx         -- not yet registered
  elseif type(ref) == 'number' then col[key]    = {ref, idx}  -- already one index
  else                              ref[#ref+1] = idx         -- already many indexes
  end
end

-- add a new row

local function add_row(self, row)
  local nrow = #self[1]+1

  if #row > 0 then -- vector of index-value
    if #row ~= #self then invalid_len(self, nrow) end
    for i,v in ipairs(row) do rawset(self[i], nrow, v) end

  else              -- vector of key-value
    local ncol, col = 0
    for k,v in pairs(row) do
      col = rawget(self, k)
      if col == nil then invalid_key(self, nrow) end
      rawset(col, nrow, v)
      ncol = ncol + 1
    end
    if ncol ~= #self then
      error("invalid number of key-value parameters in row add")
    end
  end

  if self._refcol then add_refkey(self, nrow) end
end

-- set row values

local function set_row(self, idx, val)
  local n = #val

  if n ~= 0 then    -- vector of index-value
    if n > #self then invalid_idx(self, self._refcol[idx]) end
    for i,v in ipairs(val) do self[i][idx] = v end

  else              -- vector of key-value
    for k,v in pairs(val) do
      local col = rawget(self, k)
      if col == nil then invalid_col(self, self._refcol[idx]) end
      col[idx] = v
    end
  end
end

local function get_refkey_mt(self)
  return {
    __len = function(row)
      local idx = row._idx
      return type(idx) == 'table' and #idx or 0
    end,

    __index = function(row, col)
      local idx = row._idx

      -- print('acc.__index: idx=', idx, ', key=', col)

      if type(idx) == 'number' and type(col) == 'string' then
        local vec = rawget(self, col) or invalid_col(self, self._refcol[idx])
        return      rawget(vec,  idx) or invalid_row(self, self._refcol[idx])
      end

      if type(idx) == 'table' and type(col) == 'number' then
        row._idx = idx[col] or invalid_cnt(self, self._refcol[idx[1]])
        return row
      end

      invalid_cnt(self, self._refcol[type(idx) == 'table' and idx[1] or idx])
    end,

    __newindex = function(row, col, val)
      local idx = row._idx

      -- print('acc.__newindex: idx=', idx, ', key=', col, ', val=', val)

      if type(idx) == 'number' and type(col) == 'string' then
        local vec = rawget(self, col) or invalid_col(self, self._refcol[idx])
        if type(val) == 'table' then     invalid_val(self, self._refcol[idx]) end
        rawset(vec, idx, val)
        return
      end

      if type(idx) == 'table' and type(col) == 'number' then
        idx = idx[col] or            invalid_cnt(self, self._refcol[idx[1]])
        if type(val) ~= 'table' then invalid_tbl(self, self._refcol[idx   ]) end
        set_row(self, idx, val)
        return
      end

      invalid_cnt(self, self._refcol[type(idx) == 'table' and idx[1] or idx])
    end
  }
end

local function clr_refcol(self)
  self._refcol, self._refidx = nil
  self._refkey, self._refkey_mt = {}
end

local function set_refcol(self, rcol, ridx)
  self._refcol, self._refidx = rcol, ridx
  self._refkey, self._refkey_mt = {}, get_refkey_mt(self)
  for i=1,#rcol do add_refkey(self, i) end
end

-- initialization

local function init(cols, name)
  local self = { _header = {}, _colnames = {} }
  local rcol, ridx = nil

  -- create header
  M.set_key(self, { name=name })
  M.set_key(self, { type='no-type' })
  M.set_key(self, { title='no-title' })
  M.set_key(self, { origin='no-origin' })
  M.set_key(self, { date=os.date"%d/%m/%y" })
  M.set_key(self, { time=os.date"%H:%M:%S" })

  -- create columns
  for i,v in ipairs(cols) do
    self[i] = {}
    if is_list(v) then rcol = self[i] ; ridx = i ; v = v[1] end
    self[v] = self[i]
    self._colnames[i] = v
  end

  if rcol then set_refcol(self, rcol, ridx) else clr_refcol(self) end
  return self
end

-- methods -------------------------------------------------------------------o

function M:get_key(key)
  return self._header[key]
end

function M:set_key(keys)
  local hdr = self._header
  for k,v in pairs(keys) do
    if not hdr[k] then hdr[#hdr+1] = k end
    hdr[k] = v
  end
  return self
end

function M:get_length()
  return #self[1]
end

function M:get_column_names()
  return self._colnames
end

function M:clr_refcol()
  clr_refcol(self)
  return self
end

function M:set_refcol(rcol)
  if not (type(rcol) == 'string' and self[rcol]) then
    error('invalid reference column name')
  end
  set_refcol(self, self[rcol])
  return self
end

function M:add_col(name, values)
  -- todo
  return self
end

function M:add_row(row)
  add_row(self, row)
  return self
end

function M:write(filename, columns)

  -- open file
  local name = filename or self:get_key'name' or 'tmptable'
  local file, err = io.open(name .. '.tfs', 'w')
  if not file then
    error("unable to open file '"..name.."' for writing: "..err) 
  end

  -- dump header
  local hdr = self._header
  for i=1,#hdr do
    local k, v = hdr[i], hdr[hdr[i]]
    file:write(string.format('@ %-18s %%%02ds "%s"\n', k, #v, v))
  end

  -- dump col names
  local cols = columns or self._colnames
  local ncol, nrow = #cols, #self[cols[1]]

  file:write('*')
  for icol=1,ncol do
    file:write(string.format(' %-17s ', cols[icol]))
  end
  file:write('\n')

  -- dump col types
  file:write('$')
  for icol=1,ncol do
    local fmt = type(self[cols[icol]][1]) == "number" and '%le' or '%s'
    file:write(string.format(' %-17s ', fmt))
  end
  file:write('\n')

  -- dump rows
  for irow=1,nrow do
  file:write(' ')
  for icol=1,ncol do
    local e = self[cols[icol]][irow]
    if type(e) == "number" then
      file:write(string.format('% -18.10g ', e))
    else
      file:write(string.format('%-18s ', '"'..e..'"'))
    end
  end
  file:write('\n')
  end

  -- close file
  file:close()
  return self
end

-- metamethods ---------------------------------------------------------------o

-- constructor of lines, can be unamed (inherit its name)
function MT:__call(a)
  if type(a) == 'string' then
    return function(t)
      if is_list(t) then
        return setmetatable(init(t, a), self)
      end
      error ("invalid table constructor argument, list expected")
    end
  end

  if is_list(a) then
    return setmetatable(init(a), self)
  end
  error ("invalid table constructor argument, string or list expected")
end

function M.__add(self, a)
  return self:add_row(a)
end

function M.__index(self, key)
  local idx = self._refkey[key]
  if idx == nil then
    return getmetatable(self)[key] or invalid_row(self, key)
  end

  -- print('tbl.__index: idx=', idx, ', key=', key)

  return setmetatable({_idx=idx}, self._refkey_mt)
end

function M.__newindex(self, key, val)
  local idx = self._refkey[key]

  -- print('tbl.__newindex: idx=', idx, ', key=', key, ', val=', val)

  if      idx  ==  nil    then invalid_row(self, key) end
  if type(idx) == 'table' then invalid_cnt(self, key) end
  if type(val) ~= 'table' then invalid_tbl(self, key) end
  set_row(self, idx, val)
end

-- end -----------------------------------------------------------------------o
setmetatable(_G, _G_mt) ; _G_mt = nil

return M
