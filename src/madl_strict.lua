-- The detecting of undeclared vars is discussed on:
-- http://www.lua.org/pil/14.2.html
-- http://lua-users.org/wiki/DetectingUndefinedVariables

local IGNORED_EXTRAS
local IGNORED_WRITES = {}
local IGNORED_READS  = {
  _PROMPT =true,
  _PROMPT2=true,
}

local MT = {
  __index = function (table, key)
    if IGNORED_READS[key] then return end
    local info = debug.getinfo(2, "Sl")
    MAD.warn("%s:%s: attempt to read undeclared global variable: %s\n",
             tostring(info.short_src), tostring(info.currentline), key)
    return nil
  end,

  __newindex = function (table, key, value)
    if not (IGNORED_WRITES[key] or IGNORED_EXTRAS[key]) then
      local info = debug.getinfo(2, "Sl")
      MAD.warn("%s:%s: write to undeclared global variable: %s\n",
               tostring(info.short_src), tostring(info.currentline), key)
    end
    rawset(table, key, value)
  end,
}

local require_orig = require

local function require_strict (modname)
  IGNORED_WRITES[modname] = true
  return require_orig(modname)
end

-- Raises an error when an undeclared variable is read or written.
local STRICTMODE = false

local function strict_unset ()
  local mt = getmetatable(_G)
  if mt ~= nil and mt ~= MT then
    error("invalid global metatable (i.e. not set by strict)")
  end
  setmetatable(_G, nil)
  IGNORED_EXTRAS, STRICTMODE = nil, false
  require = require_orig
end

local function strict_set (extra)
  local mt = getmetatable(_G)
  if mt ~= nil and mt ~= MT then
    error("a global metatable already exists")
  end
  if extra == nil or extra == true then extra = {} end
  assert(type(extra) == "table", "invalid argument #1 (table expected)")
  setmetatable(_G, MT)
  IGNORED_EXTRAS, STRICTMODE = extra, true
  require = require_strict
end

local function strict (status)
  local mode = STRICTMODE

      if status == false    then strict_unset()
  elseif status ~= 'status' then strict_set(status) end

  return mode
end

return { strict = strict }
