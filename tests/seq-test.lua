sequ = require 'sequence'
line = require 'line'
elem = require 'element'

rbend = elem.rbend
quad  = elem.quadrupole
drift = elem.drift

--[[
mb = sbend 'mb' {}
mq = quad  'mq' {}

lhc_cell = sequ 'lhc_cell' { 2*mb, mq }

collectgarbage'collect'
m0 = collectgarbage'count'
t0=os.clock()

lhc_seq  = sequ 'lhc_seq' { lhc_cell, { _rep=10, -2*lhc_cell, 2*lhc_cell }, -lhc_cell }

print('timing = ', t0-os.clock())
collectgarbage'collect'
print('memory = ', collectgarbage'count' -m0)

sbend:show(); io.write('\n')
quad:show(); io.write('\n')
mb:show(); io.write('\n')
mq:show(); io.write('\n')
lhc_cell:show()
lhc_seq:show()
lhc_seq:remove('mq')
lhc_seq:show()
]]

--[[
QF:QUADRUPOLE,...; // focusing quadrupole
QD:QUADRUPOLE,...; // defocusing quadrupole
B1:RBEND,...; // bending magnet of type 1
B2:RBEND,...; // bending magnet of type 2
DS:DRIFT,...; // short drift space
DM:DRIFT,...; // drift space replacing two bends
DL:DRIFT,...; // long drift space
]]

qf = quad  'qf' { length=1 }
qd = quad  'qd' { length=1 }
b1 = rbend 'b1' { length=3 }
b2 = rbend 'b2' { length=3 }
ds = drift 'ds' { length=1 }
dm = drift 'dm' { length=2  , rigid=true }
dl = drift 'dl' { length=0.5, rigid=true }

--[[
// The SPS machine is represented by the lines
SPS:   LINE=(6*SUPER);
SUPER: LINE=(7*P44,INSERT,7*P44);
INSERT:LINE=(P24,2*P00,P42);
P00:   LINE=(QF,DL,QD,DL);
P24:   LINE=(QF,DM,2*B2,DS,PD);
P42:   LINE=(PF,QD,2*B2,DM,DS);
P44:   LINE=(PF,PD);
PD:    LINE=(QD,2*B2,2*B1,DS);
PF:    LINE=(QF,2*B1,2*B2,DS);
]]

--collectgarbage'collect'
--m0 = collectgarbage'count'
--t0=os.clock()

tt = qf + 2*b1 + 2*b2 + ds
pf =line 'pf'  {qf,2*b1,2*b2,ds}
pd =line 'pd'  {qd,2*b2,2*b1,ds}
p24=line 'p24' {qf,dm,2*b2,ds,pd}
p42=line 'p42' {pf,qd,2*b2,dm,ds}
p00=line 'p00' {qf,dl,qd,dl}
p44=line 'p44' {pf,pd}

insert=line 'insert' {p24,2*p00,p42}
super =line 'super'  {7*p44,insert,7*p44}
sps   =sequ 'sps'    {6*super}

--print('timing = ', os.clock()-t0)
--collectgarbage'collect'
--print('memory = ', collectgarbage'count'-m0)

sps:show{{'length','l'}, {'kind','type'}}
--[[
seq_pf = sequ 'pf' { pf }
seq_pf:show{{'length','l'}, {'kind','typ'}}
seq_tt = sequ 'tt' { tt }
seq_tt:show{{'length','l'}, {'kind','typ'}}
]]
