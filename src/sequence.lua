--[=[
 o----------------------------------------------------------------------------o
 |
 | Sequence module
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
  sequence -- build sequences

SYNOPSIS
  sequ = require "sequence"
  my_seq = sequ ['name'] { element_list... }

DESCRIPTION
  The module mad.sequence creates new sequences and lines supported by MAD.
  The elements are not copied but referenced, i.e. store the orginal ones.

RETURN VALUE
  The object (table) that represents the flat sequence.

EXAMPLE
  sequ = require "sequence"
  elem = require "element"
  MB, MQ = elem.sbend, elem.quadrupole
  my_seq = sequ 'name' {
    MQ 'QF', MB 'MB', MQ 'QD', MB 'MB',
  }

SEE ALSO
  line, element, beam
]]

M.__help.add = [=[
  :add(item [, pos [, from [, refer]]])                      :positional params
  :add{[item=]item [, at=pos [, from='from' [, refer='refer']]]}  :named params
  item= elem|list|sequ
  from= 'start'|'end'|'prev'|'next'|'name'  (default= at and 'start' or 'prev')
  refer= 'entry'|'centre'|'exit'|'name'        ('name' is only for subsequence)
]=]

-- requires ------------------------------------------------------------------o

local object = require 'object'
local utils  = require 'utils'
local line   = require 'line'
local marker = require 'element'.marker

-- locals --------------------------------------------------------------------o

local getmetatable, setmetatable = getmetatable, setmetatable
local type, ipairs = type, ipairs
local is_list, show_list = utils.is_list, utils.show_list

-- metatable for the root of all sequences
local MT = object {} 

 -- make the module the root of all sequences
MT (M)
M.name = 'sequence'
M.kind = 'sequence'
M.is_sequence = true

-- sequence fields
local sequence_fields = { 'length', 'refer' }

-- sequence markers (unique instances)
local start_marker = marker '$start' {}
local end_marker   = marker '$end'   {}

-- functions -----------------------------------------------------------------o

-- utils
local function copy_fields(t, a, lst)
  for _,v in ipairs(lst) do
    t[v] = a[v]
  end 
end

local function show_fields(t, lst)
  local a, k, s
  lst = lst or sequence_fields
  for _,v in ipairs(lst) do
    if is_list(v) then k, s = v[1], v[2] else k, s = v, v end
    a = t[k]
    if a then io.write(', ', s, '= ', tostring(a)) end
  end
end

-- search
local function find_index_by_ref(t, a, start)
  for i=start or 1,#t do -- linear search
    if t[i] == a then return i end
  end
end

local function find_index_by_name(t, name, start)
  for i=start or 1,#t do -- linear search
    if t[i].name == name then return i end
  end
  return nil
end

local function find_index_by_pos(t, s_pos, start)
  for i=start or 1,#t do  -- binary search: TODO
    if t[i].s_pos >= s_pos then return i end
  end
  return nil
end

-- dictionnary
-- shadow elements (special inheritance for shared elements)
local function shadow_get(elem)
  return getmetatable(elem).__index
end

local function shadow_class(elem)
  return shadow_get(elem):class()
end

local function shadow_element(elem, idx, at, from, refer, seq)
  local states = {s_pos='todo', _idx=idx, _at=at, _from=from, _refer=refer, _seq=seq,
                  element=shadow_get, class=shadow_class}

  if elem.shared then -- read-write element
    return setmetatable(states, {__index=elem, __newindex=elem})
  else                -- read-only element
    return elem(states)
  end
end

local function clean_shadow(self)
  for _,v in ipairs(self) do
    v._idx, v._at, v._from, v._refer, v._seq = nil
  end
end

-- construction
local function add_element_key(self, elem)
  if elem.kind == 'drift' then return end     -- drifts are never registered

-- io.write('add_elem: \'', elem.name, '\' at slot ', #self+1, '\n')
  local name = elem.name
  local ref = self[name]

      if ref == nil     then self[name] = elem          -- not yet registered
  elseif ref.is_element then self[name] = {ref, elem}   -- already one element
  else                       ref[#ref+1] = elem         -- already many elements  
  end
end

local function add_element(self, elem, at, from, refer, seq)
  local i = #self+1
  self[i] = shadow_element(elem, i, at, from, refer, seq)
  add_element_key(self, self[i])
end

local function add_sequence(self, sequ, at, from, refer, rev)
--  io.write('add_sequ: \'', sequ.name, '\' at slot ', #self+1, '\n')
  if rev and rev<0 then -- reverse, store sequ info with end_marker
    add_element(self, shadow_get(sequ[#sequ]), at, from, refer, sequ)
    for j=#sequ-1,1,-1 do 
      at = sequ[j+1].s_pos - sequ[j].s_pos - sequ[j].length
      add_element(self, shadow_get(sequ[j]), at, 'prev', 'entry')
    end
  else                  -- direct, store sequ info with start_marker
    add_element(self, shadow_get(sequ[1]), at, from, refer, sequ)
    for j=2,#sequ do
      at = sequ[j].s_pos - sequ[j-1].s_pos - sequ[j-1].length
      add_element(self, shadow_get(sequ[j]), at, 'prev', 'entry')
    end
  end
end

local add_line -- forward declaration (x-ref in this order lets add_item to be inlined)

local function add_error(self, item)
  error("invalid item '"..(item.name or '?').."' at slot "..(#self+1).." in sequence '"..self.name.."'")
end

local function add_item(self, item, at, from, refer, rev)
      if item.is_element  then add_element (self, item, at, from, refer)
  elseif item.is_line     then add_line    (self, item, at, from, refer, rev)
  elseif item.is_sequence then add_sequence(self, item, at, from, refer, rev)
  else                         add_error   (self, item)
  end
end

add_line = function(self, line, at, from, refer, rev)
--  io.write('add_line: \'', line:mangled_name(), '\' at slot ', #self+1, '\n')
  local j_beg, j_end, j_step
  local n = (line._rep or 1) * (rev or 1)

  if n<0
  then n, j_beg, j_end, j_step = -n, #line, 1, -1
  else    j_beg, j_end, j_step =     1, #line,  1
  end

  for i=1,n do
    for j=j_beg,j_end,j_step do
      add_item(self, line[j], at, from, refer, j_step)
    end
  end
end

local function add_list(self, lst)
  if is_list(lst) then
    for _,v in ipairs(lst) do
      add_item(self, v)
    end
  else add_item(self, lst)
  end
end

-- compute s positions of sequence elements
local function spos_error(self, elem)
  if elem.s_pos == 'ongoing' then
    error('cycling dependencies detected in sequence '..self.name..' for element '..elem.name)
  else
    error('invalid element detected in sequence '..self.name..' for element '..elem.name)
  end
end

local function element_spos(self, elem)
  local s_pos = elem.s_pos
  if type(s_pos) == 'number' then return s_pos end
  if s_pos ~= 'todo' then spos_error(self, elem) end

  elem.s_pos = 'ongoing'

  local len   = self.length
  local idx   = elem._idx
  local seq   = elem._seq
  local pos   = elem._at    or 0
  local from  = elem._from  or elem._at and 'start' or 'prev'
  local refer = elem._refer or self.refer           or 'entry'

      if refer == 'entry'   then ;
  elseif refer == 'centre'  then pos = pos - (seq and seq.length/2 or elem.length/2)
  elseif refer == 'exit'    then pos = pos - (seq and seq.length   or elem.length)
  elseif seq and seq[refer] then pos = pos - seq[refer].s_pos
  else error("invalid refer: "..elem._refer)
  end

      if from == 'start'       then ;
  elseif from == 'end' and len then pos = len - pos
  elseif from == 'prev'        then pos = pos + self[idx-1].s_pos + self[idx-1].length
  elseif from == 'next'        then pos = pos + element_spos(self, self[idx+1])
  elseif self[from]            then pos = pos + element_spos(self, self[from])
  else error("invalid from: "..elem._from)
  end

  elem.s_pos = pos
  return pos
end

local function sequence_spos(self)
  for _,v in ipairs(self) do
    element_spos(self, v)
  end
end

local function adjust_length(self)
  local em  = self[#self] -- end_marker
  local len = self.length
  if not len or len < em.s_pos then self.length = em.s_pos
  elseif        len > em.s_pos then em.s_pos = len
  end
end

local function clean_sequence(self)
  local j=2 -- keep start_marker
  for i,v in ipairs(self) do
    if shadow_get(v) ~= start_marker and
       shadow_get(v) ~= end_marker   and
      (v.kind ~= 'drift' or v.rigid) then self[j], j = self[i], j+1 end
  end
  self[j] = self[#self] -- copy end_marker
  for i=j+1,#self do self[i] = nil end -- clean remaining slots

  clean_shadow  (self)
end

-- methods -------------------------------------------------------------------o

function M:show(disp)
  io.write("sequence '", self.name,"' { ")
  show_list(self, disp)
  io.write(' }\n')
  for _,v in ipairs(self) do v:show(disp) end
  io.write('endsequence\n')
end

function M:show_madx(disp)
  io.write("'", self.name, "': sequence, ")
  show_list(self, disp)
  io.write(';\n')
  for _,v in ipairs(self) do v:show_madx(disp) end
  io.write('endsequence;\n')
end

function M:add(a, at, from, refer)
  if at or not a.at then add_item(self, a, at, from, refer)                 -- positional params
  elseif is_list(a) then add_item(self, a.item or a[1], a.at, a.from, a.refer)   -- named params
  else error("invalid set of parameters in incremental construction of sequence '"..self.name.."'")
  end
  return self
end

function M:set(a)
  self[1] = nil -- clear sequence
  copy_fields(self, a, sequence_fields)
  add_element(self, start_marker, 0)
  add_list   (self, a)
  return self:done()
end

function M:done()
  add_element   (self, end_marker, self.length)
  sequence_spos (self)
  adjust_length (self)
  clean_sequence(self)
  return self
end

-- metamethods ---------------------------------------------------------------o

-- constructor of sequences, can be unamed (inherit its name)
function MT:__call(a)
  if type(a) == 'string' then
    return function(t)
      if is_list(t) then
        self.__index = self         -- inheritance
        return setmetatable({name=a}, self):set(t)
      end
      error ("invalid sequence constructor argument, list expected")
    end
  end

  if is_list(a) then
    self.__index = self             -- inheritance
    return setmetatable({}, self):set(a)
  end
  error ("invalid sequence constructor argument, string or list expected")
end

-- construction
function M.__add(sequ, a)
  return sequ:add(a)
end

-- repetition
function M.__mul(n, sequ)
  if type(sequ) == 'number' then n, sequ = sequ, n end
  return line { _rep=n, sequ }
end

-- reflection
function M.__unm(sequ, _)
  return line { _rep=-1, sequ }
end 

-- end -------------------------------------------------------------------------
return M


--[=[
--------------------------------------------------------------------------------
-- TODO ------------------------------------------------------------------------
--------------------------------------------------------------------------------

--[[
local function flatten(self, name)
  local t = { name=name or self.name }

  for i,v in ipairs(self) do
    if type(v) == 'function' then v = v() end
    if v.is_sequence then
      add_sequence(t, v)
    elseif v.is_element then
      add_element(t, v)
    elseif v.is_line then
      add_line(t, v)
    else
      error('invalid sequence element at slot '..i)
    end
  end

  return localise(t)
end
--]]

-- geometry -- TODO: at, from, refer
local function localise(self, start)
  local s_pos = start and v[start].s_pos or 0
  for i=start or 1,#self do
    self[i].s_pos = s_pos
    if self[i].length then s_pos = s_pos + self[i].length end
  end
  if not self.length then self.length = s_pos end
  return self
end

local function insert_element_key(self, elem)
  local name = elem.name              -- dict part
  local ref = self[name]
  if ref == nil then                  -- not yet registered
    self[name] = elem
  elseif ref.is_element then          -- already one element
    self[name] = ref.i_pos < elem.i_pos and {ref, elem} or {elem, ref}
  else                                -- already many elements
    table.insert(ref, find_index_by_idx(ref, elem.i_pos), elem)
  end
end

local function remove_element_key(self, elem)
  local name = elem.name              -- dict part
  local ref = self[name]
  if ref.is_element then              -- single element
    self[name] = nil
  else                                -- list of elements
    table.remove(ref, find_index_by_ref(ref, elem))
    if #ref == 1 then self[name] = ref[1] end -- was a pair
  end
end

-- edition -- TODO: check s_pos and length
local function insert_element(self, elem, before)
  test_membership(self, before)
  local i = before.i_pos
  table.insert(self, i, elem)
  update_index(self, i)
  insert_element_key(self, elem)
end

local function remove_element(self, elem)
  test_membership(self, elem)
  local i = elem.i_pos
  remove_element_key(self, elem)
  table.remove(self, i)
  update_index(self, i)
end

local function replace_element(self, old_elem, new_elem)
  test_membership(self, old_elem)
  local i = old_elem.i_pos
  self[i] = new_elem
  new_elem.i_pos = i
  remove_element_key(self, elem)
  insert_element_key(self, elem)
end

local function swap_elements(self, elem1, elem2, update_key)
  test_membership(self, elem1)
  test_membership(self, elem2)
  local i1, i2 = elem1.i_pos, elem2.i_pos
  self[i1], self[i2] = elem2, elem1
  elem1.i_pos, elem2.i_pos = i2, i1
  if update_key then
    remove_element_key(self, elem1)
    remove_element_key(self, elem2)
    insert_element_key(self, elem1)
    insert_element_key(self, elem2)
  end
end

function M:remove(a, count) -- TODO
  if type(a) == 'string' then
    a = self[a]
    if is_list(a) then a = a[count or 1] end
  end
  remove_element(self, a)
end

function M:insert(a, at, count) -- TODO
  if type(at) == 'string' then
    at = self[at]
    if is_list(at) then at = at[count or 1] end
  elseif type(at) == 'number' then
    at = self[find_index_by_pos(self, at)]
  end
  insert_element(self, a, at)
end
--]=]
