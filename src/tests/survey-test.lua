-- build sequence

sequ = require 'sequence'
elem = require 'element'

mid  = elem.marker 'midseq' {}
dip  = elem.sbend 'dip'     { length=10, angle=2*math.pi/12 }
dipt = dip 'dipt'           { tilt=0.1 }

musr = sequ 'musr' { refer=centre, length=120 }

musr
:add { dip'dip1' {}, at=5}
:add ( dip'dip2' {})  
:add ( dip'dip3' {})  
:add ( dip'dip4' {})  
:add ( dip'dip5' {})  
:add ( dip'dip6' {})  

:add ( mid, 60)  

:add ( dip 'dip7'  {})  
:add ( dipt'dip8'  {})  
:add ( dipt'dip9'  {})  
:add ( dipt'dip10' {})  
:add ( dipt'dip11' {})
:add ( dip 'dip12' {})
:done()

-- compute survey

survey=require 'survey'.survey

io.write("computing survey\n")
tbl = survey { seq=musr, tbl='survey' }

io.write("writing survey table\n")
tbl:write()