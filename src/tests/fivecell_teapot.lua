local fivecell = require 'fivecell' .thin()

-- output madx-like sequence
fivecell:show_madx{{'length','l'}}

-- compute survey
local survey = require 'survey' .survey

io.write("computing survey\n")
local tbl = survey { seq = fivecell, tbl = 'survey' }

io.write("writing survey table\n")
tbl:write()

-- build the tacked map

local track = require 'track' .track
local map   = require 'mflow' .map

if false then
  m = map {v={'x','px', 'y','py', 't','pt'}, mo={0,0,0,0,0,0}} -- only scalars

  m.x  = 1e-4 -- orbit at origin
  m.y  = 1e-4
  m.t  = 0

  m.px = 0
  m.py = 0
  m.pt = 1e-6
else
  m = map {v={'x','px', 'y','py', 't','pt'}, mo={2,2,2,2,2,2}} -- high orders

m:print()

  m.x:set({0}, 1e-4)
  m.y:set({0}, 1e-4)
--  m.t:set({0}, 0   )

  m.px:set({0}, 0)
  m.py:set({0}, 0)
--  m.pt:set({0}, 1e-6)   -- delta p

  m.x :set({1,0,0,0,0,0}, 1) -- derivatives
  m.px:set({0,1,0,0,0,0}, 1)
  m.y :set({0,0,1,0,0,0}, 1)
  m.py:set({0,0,0,1,0,0}, 1)
--  m.t :set({0,0,0,0,1,0}, 1)
--  m.pt:set({0,0,0,0,0,1}, 1)

-- scalars case
  m.t  = 0
  m.pt = 1e-6

end

m:print()
m.charge, m.dir, m.total_path, m.beta0_inv = 1, 1, 0, 1.0000001305599685659
tbl = track { seq = fivecell, map = m, tbl = 'track' }
m:print()
tbl:write()
