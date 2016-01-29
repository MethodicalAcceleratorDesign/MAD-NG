local function make_fodo()
  local sequ = require 'sequence'
  local elem = require 'element'
  local line = require 'line'
  local drift = elem.drift
  local sbend = elem.sbend
  local quad  = elem.quadrupole
  local sext  = elem.sextupole
  local pi = math.pi
  local d1,d2,bd,qf,qd,sf,sd

  qf = quad'qf' { length=0.5, k1=  0.6536 }
  qd = quad'qd' { length=0.5, k1= -0.6556 }

  bd = sbend'bd' { length=1.5, angle= 2*pi/100 }

  sf = sext'sf' { length=0.1, k2=  24.9056 }
  sd = sext'sd' { length=0.1, k2= -48.4914 }

  d1 = drift'd1' { length= .25, rigid=true } 
  d2 = drift'd2' { length=0.15, rigid=true }

  cell = line'cell' { qf,sf,d2,bd,d1,qd,qd,sd,d2,bd,d1,qf } 

  return sequ 'madsync' { 50*cell }
end

madsync = make_fodo()
madsync:show_madx{ {'length','l'}, 'angle', 'k1', 'k2' }
