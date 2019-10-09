--[=[
 o-----------------------------------------------------------------------------o
 |
 | Lua Hook function (from Lua repository)
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
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
  - Provide a standard interface to trace Lua code.

 o-----------------------------------------------------------------------------o
]=]

local function hook(event, newlines_)
  -- event: "call", "tail call", "return", "line", "count"
  local t = debug.getinfo(3)

  io.write(" >>> ")
  if t ~= nil and t.currentline >= 0 then
    local src = t.short_src:match("[^/]*$")
    io.write(src, ":", t.currentline, " ")
  end

  t = assert(debug.getinfo(2), "unexpected nil context")
  local src = t.short_src:match("[^/]*$")

  if t.what == "main" then
    if event == "call" then
      io.write("begin ", src)
    elseif event == "return" then
      io.write("end "  , src)
    end
  elseif t.what == "Lua" then
    io.write(event, " ", t.name or "(Lua)", " <", t.linedefined, ":", src, ">")
  elseif t.what == "C" then
    io.write(event, " ", t.name or "(C)", " [", t.what, "] ")
  else
    io.write("unknown what: ", t.what, " for event: ", event)
  end

  io.write("\n")
end

local function dbghook (mode_, count_)
  -- mode_: "c", "r", "l" or "off"
  if mode_ == "off"
  then debug.sethook() -- turn hook off
  else debug.sethook(hook, mode_ or "cr", count_)
  end
end

return { dbghook = dbghook }

