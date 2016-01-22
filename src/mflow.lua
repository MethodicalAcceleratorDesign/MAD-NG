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

local gmath   = require 'gmath'
local tpsa    = require 'tpsa'
local xtpsa   = require 'xtpsa'
local complex = require 'complex'

local desc  = xtpsa.desc
local mono  = xtpsa.mono

-- locals --------------------------------------------------------------------o

local type, insert = type, table.insert
local isnum, iscpx, isscl, ident =
      gmath.is_number, gmath.is_complex, gmath.is_scalar, gmath.ident

local D = {}  -- private desc
local V = {}  -- private keys
local T = {}  -- temporary keys

local S_lst = {'x', 'px', 'y', 'py', 't', 'pt'} -- allowed variable names
local S_dft = {'x', 'px', 'y', 'py', 't', 'pt'} -- default variable names

-- implementation ------------------------------------------------------------o

local function make_map(args, num_ctor, tpsa_ctor)
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
      map[V][var[i]] = num_ctor(0)
    else
      map[V][var[i]] = tpsa_ctor(dsc, mo[i])
    end
  end

  return setmetatable(map, M)
end

local function map(args)
  return make_map(args, ident, tpsa)
end

local function cmap(args)
  return make_map(args, complex, ctpsa)
end

-- indexing ------------------------------------------------------------------o

function M.__index (tbl, key)
--  io.write("getting ", key, '\n')
  return tbl[V][key] or tbl[T][key] or M[key]
end

function M.__newindex (tbl, key, val)
--  io.write("setting ", key, '\n')
  local K = tbl[V][key] and V or T
  local v = tbl[K][key]         -- get the variable

  if v == nil then              -- create the variable
    tbl[K][#tbl[K]+1] = key
    tbl[K][key]       = val
  elseif isscl(v) then          -- number or tpsa -> number
    tbl[K][key] = isscl(val) and val or val.coef[0]
  elseif isscl(val) then
    v:scalar(val)               -- number -> TPSA
  else
    val:copy(v)                 -- TPSA -> TPSA
  end
end

-- TODO

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

function M.to(tbl, ...)
  if not tbl[D] then return end
  local to, mo = -1
  for _,v in ipairs{...} do
    mo = isnum(v) and 0 or v.mo
    to = mo > to and mo or to
  end

  tpsa.gtrunc(tbl[D],to)
end

function M.set(tbl, var, mono, val)
  if isnum(tbl[var]) then
    assert(mono_sum(mono) == 0, "Invalid set for constant var")
    tbl[var] = val
  else
    tbl[var]:set(mono, val)
  end
end

function M.get(tbl, var, mono)
  return isnum(tbl[var]) and assert(mono_sum(mono) == 0, "Invalid get") and tbl[var]
         or tbl[var]:get(mono)
end

------------------------------------------------------------------------------o

return {
  map  = map,
  cmap = cmap
}
