sequ = require "sequence"
elem = require "element"
line = require "line"

local function make_als()
  local L1,L2,L3,L4,L5,L6,L7,L8,L9,L10
  local L11,L12,L13,L14,L15,L16,L17,L18,L19,L20
  local L21,L22,L23,L24,L25,L26,L27,L27A,L27B,L27C,L27D,DS
  local QF1,QF2,QD1,QD2,QFA1,QFA2,SF,SD,CAV,BEND,VC5,BEND1
  local SUP1,SUPB

  local drift = elem.drift
  local rbend = elem.rbend
  local quad  = elem.quadrupole
  local sext  = elem.sextupole
  local mark  = elem.marker
  local rfcav = elem.rfcavity

	L1   = drift "L1"     { length=2.832695 }
  L2   = drift "L2"     { length=0.45698 }
  L3   = drift "L3"     { length=0.08902 }
  L4   = drift "L4"     { length=0.2155 }
  L5   = drift "L5"     { length=0.219 }
  L6   = drift "L6"     { length=0.107078 }
  L7   = drift "L7"     { length=0.105716 }
  L8   = drift "L8"     { length=0.135904 }
  L9   = drift "L9"     { length=0.2156993 }
  L10  = drift "L10"    { length=0.089084 }
  L11  = drift "L11"    { length=0.235416 }
  L12  = drift "L12"    { length=0.1245 }
  L13  = drift "L13"    { length=0.511844 }
  L14  = drift "L14"    { length=0.1788541 }
  L15  = drift "L15"    { length=0.1788483 }
  L16  = drift "L16"    { length=0.511849 }
  L17  = drift "L17"    { length=0.1245 }
  L18  = drift "L18"    { length=0.235405 }
  L19  = drift "L19"    { length=0.089095 }
  L20  = drift "L20"    { length=0.2157007 }
  L21  = drift "L21"    { length=0.177716 }
  L22  = drift "L22"    { length=0.170981 }
  L23  = drift "L23"    { length=0.218997 }
  L24  = drift "L24"    { length=0.215503 }
  L25  = drift "L25"    { length=0.0890187 }
  L26  = drift "L26"    { length=0.45698 }
  L27  = drift "L27"    { length=2.832696 }
  L27a = drift "L27a"   { length=0.8596 }
  L27b = drift "L27b"   { length=0.1524 }
  L27c = drift "L27c"   { length=0.04445 }
  L27d = drift "L27d"   { length=1.776246 }
  DS   = drift "DS"     { length=0.1015 }

  QF1  = quad "QF1"     { length=0.344, k1= 2.2474+6.447435260914397e-03 }
  QF2  = quad "QF2"     { length=0.344, k1= 2.2474 }
  QD1  = quad "QD1"     { length=0.187, k1=-2.3368-2.593018157427161e-02 } 
  QD2  = quad "QD2"     { length=0.187, k1=-2.3368 }  
  QFA1 = quad "QFA1"    { length=0.448, k1= 2.8856 }  
  QFA2 = quad "QFA2"    { length=0.448, k1= 2.8856 }  

  local ksf= -41.67478927130080+0.3392376315938252
  local ksd=  56.36083889436033-0.1043679358857811

  SF= sext "SF" { length=2*0.1015, k2=ksf }
  SD= sext "SD" { length=2*0.1015, k2=ksd }

  VC5 = mark "VC5" {}

  local l_bend = 0.86621
  local alpha  = 0.17453292519943295769236907684886

  BEND  = rbend "BEND"  { length=l_bend, angle=alpha, k1=-0.778741 }
  BEND1 = rbend "BEND1" { length=l_bend, angle=alpha, k1=-0.778741 }

  CAV = rfcav "CAV" { length=0.0, volt=-1.0, rev_freq=500.0e6 }

  SUP1 =  L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+SF+L10+
          L11+QFA1+L12+SD+L13+
          L14+BEND+L15+L16+SD+L17+
          QFA2+L18+L19+SF+L20+BEND+L21+
          L22+QD2+L23+L24+QF2+L25+
          L26+VC5+L27

  SUPB =  L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+SF+L10+
          L11+QFA1+L12+SD+L13+
          L14+BEND+L15+L16+SD+L17+ 
          QFA2+L18+L19+SF+L20+BEND1+L21+
          L22+QD2+L23+L24+QF2+L25+
          L26+VC5+L27

  return sequ 'ALS' { 11*SUP1+SUPB+CAV }
end

als = make_als()
als:show_madx{{'length','l'}, 'k1', {'kind','typ'}}
