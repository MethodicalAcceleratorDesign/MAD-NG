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

local level = 0

local function hook_(event)
 local t = debug.getinfo(3)

 io.write(level, " >>> ", string.rep(" ",level))
 if t ~= nil and t.currentline >= 0 then
  io.write(t.short_src,":",t.currentline," ")
 end

 t = debug.getinfo(2)
 if event == "call" then
  level = level+1
 else
  level = level-1
  if level < 0 then level = 0 end
 end

 if t.what == "main" then
  if event == "call"
  then io.write("begin ", t.short_src)
  else io.write("end "  , t.short_src)
  end
 elseif t.what == "Lua" then
  io.write(event," ",t.name or "(Lua)"," <",t.linedefined,":",t.short_src,">")
 else
  io.write(event," ",t.name or "(C)"," [",t.what,"] ")
 end
 io.write("\n")
end

local function hook (mode)
  debug.sethook(hook_, mode or "cr")
  level = 0
end

return { hook = hook }

