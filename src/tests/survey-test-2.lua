-- build sequence

sequ = require 'sequence'
elem = require 'element'
surv = require 'survey'.survey -- todo call semantic

local make_sequ = function (N)
	local dip = elem.sbend'dip' { length=10, angle=2*math.pi/N, tilt=2*math.pi/N }
	return sequ'seq' { length=N*10, N*dip }
end

io.write("computing sequence\n")
musr=make_sequ(10000)

io.write("computing survey\n")
tbl = surv { seq=musr, tbl='survey-2' }

io.write("writing survey table\n")
tbl:write()

-- print(tbl.s_pos.dip[150])