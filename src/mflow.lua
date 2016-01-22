--[=[
 o----------------------------------------------------------------------------o
 |
 | Map Flow module
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
  - provides constructors and functions for map flow

 o----------------------------------------------------------------------------o
]=]

local M = { __help = {}, __test = {} }

-- help ----------------------------------------------------------------------o

M.__help.self = [[
NAME
  mflow -- Map Flow module

SYNOPSIS
  local mflow = require 'mflow'

DESCRIPTION
  The module mflow provides consistent definitions and constructors of map
  flow handling vector or reals and real GTPSAs.

RETURN VALUES
  The constructors of map flow.

SEE ALSO
  gmath, complex, matrix, cmatrix, tpsa, ctpsa
]]

-- modules -------------------------------------------------------------------o

local gmath = require 'gmath'
local tpsa  = require 'tpsa'
local xtpsa = require 'xtpsa'

local desc  = xtpsa.desc
local mono  = xtpsa.mono

-- locals --------------------------------------------------------------------o

local type, insert = type, table.insert
local isnum = gmath.is_number

local D = {}  -- private desc
local V = {}  -- private keys
local T = {}  -- temporary keys

local S_lst = {'x', 'px', 'y', 'py', 't', 'pt'} -- allowed variable names
local S_dft = {'x', 'px', 'y', 'py', 't', 'pt'} -- default variable names

-- implementation ------------------------------------------------------------o

local function make_map(args, ctor)
  if args.v then
    for _,v in ipairs(args.v) do
      assert(S_lst[v], "invalid variable name")
    end
  else args.v = S_dft end

  local dsc = desc(args)
  local map = { [D]=dsc, [V]={}, [T]={} }
  local var = args.v
  local mo  = args.mo or args.vo

  for i=1,dsc.nmv do
    map[V][i] = var[i]
    if mo[i] == 0 then
      map[V][var[i]] = 0
    else
      map[V][var[i]] = ctor(dsc, mo[i])
    end
  end

  return setmetatable(map, M)
end

local function mflow(args)
  return make_map(args, tpsa)
end

local function cmflow(args)
  return make_map(args, ctpsa)
end

M.mflow  = mflow
M.cmflow = cmflow

-- indexing ------------------------------------------------------------------o

function M.__index (tbl, key)
--  io.write("getting ", key, '\n')
  return tbl[V][key] or tbl[T][key] or M[key]
end

function M.__newindex (tbl, key, val)
--  io.write("setting ", key, '\n')
  local K   = tbl[V][key] and V or T
  local var = tbl[K][key]

  if var == nil then
    tbl[K][#tbl[K]+1] = key           -- save the name
    if isnum(val) then
      tbl[K][key] = val               -- create number
    else
      tbl[K][key] = val:set_var()     -- create TPSA
    end

  elseif isnum(var) then
    if isnum(val) then
      tbl[K][key] = val               -- number -> number
    else
      tbl[K][key] = val.coef[0]       -- TPSA -> number
      val:release()
    end
  elseif isnum(val) then
    tpsa.scalar(var, val)             -- number -> TPSA
  else
    tpsa.copy(val, var)               -- TPSA -> TPSA
    val:release()
  end
end

function M.clear(tbl)
  tbl = tbl[T]
  for i,k in ipairs(tbl) do
    tbl[k]:set_tmp():release()
    tbk[i] = nil
  end
end

function M.print_tmp(tbl)
  tbl = tbl[T]
  for i,k in ipairs(tbl) do
    print(i, k, tbl[k])
  end
end

function M.print(tbl)
  for _,name in ipairs(tbl[V]) do
    local var = tbl[V][name]
    io.write(name, ': ')
    if isnum(var) then
      print(var)
    else
      var:print()
    end
  end
end

function M:to(...)
  if not self[D] then return end
  local to, mo = -1
  for _,v in ipairs{...} do
    mo = isnum(v) and 0 or v.mo
    to = mo > to and mo or to
  end

  tpsa.gtrunc(self[D],to)
end

function M.set(m, var, mono, val)
  if isnum(m[var]) then
    assert(mono_sum(mono) == 0, "Invalid set for constant var")
    m[var] = val
  else
    m[var]:set(mono, val)
  end
end

function M.get(m, var, mono)
  return isnum(m[var]) and assert(mono_sum(mono) == 0, "Invalid get") and m[var]
         or m[var]:get(mono)
end

------------------------------------------------------------------------------o

return mflow
