local table = require 'tfstable'

-- name is the reference column for tab.row.col access
local tab = table 'survey' { {'name'}, 'x', 'y', 'z', 'theta', 'phi', 'psi' }

local tab = tab
 +{ 'drift', 0.1, 0.2, 0.5, 0, 0, 0 } -- unamed array
 +{ name='mq', x=0.2, y=0.4, z=1, theta=0, phi=0, psi=0 }
 +{ name='mb', x=0.3, y=0.5, z=1, theta=0, phi=0, psi=0 }
 +{ name='mq', x=0.2, y=0.4, z=1, theta=0, phi=0, psi=0 }

---[[
collectgarbage'collect'
m0 = collectgarbage'count'
do
	t0=os.clock()
	for i=1,1000000 do
		tab:add_row{ 'mq', i*0.5, i*0.6, 1, 0, 0, 0 }
	end
	print('timing = ', os.clock()-t0)
end
collectgarbage'collect'
print('memory = ', collectgarbage'count'-m0)
--]]

print('#tab= ', #tab[1])    -- x of 'mb'
--tab:write()         -- equivalent to tab:write"name.tfs"

-- indexed access
print('tab[2][3]= '		, tab[2][3])    -- x of 'mb'
print('tab.x[2]= '		, tab.x[2])     -- x of 'mq'
print('tab.name[2]= '	, tab.name[2])  -- 'mq'

-- READINGS --
print('READINGS')

-- non-existing element
-- print('tab.qx= '		, tab.qx)

-- single element
print('tab.mb= '		, tab.mb)
print('#tab.mb= '		, #tab.mb)
-- print('tab.mb[1].x= ', tab.mb[1].x) -- count error, ok
print('tab.mb.x= '		, tab.mb.x)
-- print('tab.mb.t= '		, tab.mb.t) -- invalid column, ok

-- duplicated element
print('tab.mq= '		, tab.mq)
print('#tab.mq= '		,#tab.mq)
-- print('tab.mq.x= '		, tab.mq.x) -- count error, ok
print('tab.mq[1].x= '	, tab.mq[1].x)
-- print('tab.mq[1].t= '	, tab.mq[1].t) -- invalid column, ok

-- WRITINGS --
print('WRITINGS')

print('tab.mq[1].x = 6')
tab.mq[1].x = 6
print(tab.mq[1].x)

print('tab.mq[1] = {...}')
tab.mq[1] = { x=0.3, y=0.5, z=1, theta=0, phi=0, psi=0 }
print(tab.mq[1].x)
print(tab.mq[1].y)
print(tab.x[2])

print('tab.mb = {...}')
tab.mb = { x=0.4, y=0.6, z=1, theta=0, phi=0, psi=0 }
print(tab.mb.x)
print(tab.mb.y)
print(tab.x[3])
