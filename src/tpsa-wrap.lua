local ffi = require"ffi"

local packages = {
  ffi  = "tpsaFFI",
  mad  = "..tpsa",
  berz = "tpsaBerz",
  yang = "tpsaYang",
  mc   = "tpsaMC"
}


local M = {}

local tpsa
local curr_loaded = {}

local folder_of_this_file = (...):match("(.-)[^%.]+$")
function M.set_package(name)
  if not packages[name] then
    error("Unrecognized package title: " .. (name or "[missing title]") ..
          ". Use one of: ffi, berz, yang, mc")
  end

  if name == "mad" then
    tpsa = require("mad.tpsa")
  else
    tpsa = require(folder_of_this_file .. packages[name])
  end

  curr_loaded =  { name=name }
  setmetatable(M, { __index = tpsa })
end

function M.init(vars,no,knobs,ko)  -- same as tpsaFFI.init
  if not curr_loaded.name then
    error("Set a package first")
  elseif curr_loaded.name == "ffi" then
    curr_loaded.t = tpsa.init(vars,no,knobs,ko)

  elseif curr_loaded.name =="mad" then
    if type(vars) == "number" then
      io.stderr:write("WARNING: mad package needs different initializer. Adjusting...\n")
      local nv = #vars
      vars = {}
      for i=1,nv do vars[i] = no end
    end
    curr_loaded.t = tpsa.init(vars,vars,no)

  else
    if type(vars) == "table" or knobs or ko then
      io.stderr:write("WARNING: Currently loaded package [" .. curr_loaded.name ..
            "] does not support generalised initialisation. Adjusting...\n")
      vars = #vars
    end

    curr_loaded.t = tpsa.init(vars,no)
  end
  curr_loaded.vars = vars
  return curr_loaded.t
end

function M.new()
  return curr_loaded.t:same()
end

function M.get_map()
  local map = {}
  if curr_loaded.name == "ffi" and type(curr_loaded.vars) == "table" then
    for i=1,#curr_loaded.vars do
      map[i] = curr_loaded.t:new(curr_loaded.vars[i])
    end
  else
    for i=1,curr_loaded.vars do
      map[i] = curr_loaded.t:same()
    end
  end

  return map
end


return M
