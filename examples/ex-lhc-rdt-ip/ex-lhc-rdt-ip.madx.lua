option, warn,info;
system,"ln -s /afs/cern.ch/eng/acc-models/lhc/2022/acc-models-lhc";


beam,  bv= 1,
 particle=proton, charge=1, mass=0.938272046,
 energy= 450,   npart=1.2e11,kbunch=2556,
 ex=5.2126224777777785e-09,ey=5.2126224777777785e-09;

Option, -echo, -warn,-info;
call,file="acc-models-lhc//lhc.seq";
call,file="acc-models-lhc/operation/optics/R2022a_A11mC11mA10mL10m.madx";

Option, echo,warn,info;

use, sequence=lhcb1;


use, sequence=lhcb1;
SELECT, FLAG=makethin, CLASS=rbend, slice=0;
SELECT, FLAG=makethin, CLASS=quadrupole, slice=0;
SELECT, FLAG=makethin, CLASS=sbend, slice=0;
MAKETHIN, SEQUENCE=lhcb1;
use, sequence=lhcb1;


select, flag=twiss, pattern=BPM, column= name,s, betx,bety, mux, muy, dx, dy,x,y;
select, flag=twiss, pattern=IP, column= name,s, betx,bety, mux, muy, dx, dy,x,y;



ko=kmax_MO/Imax_MO * 40 / (450*3.33);
kof.a81b1 := ko;
kof.a12b1 := ko;
kof.a23b1 := ko;
kof.a34b1 := ko;
kof.a45b1 := ko;
kof.a56b1 := ko;
kof.a67b1 := ko;
kof.a78b1 := ko;

kod.a81b1 := ko;        
kod.a12b1 := ko;       
kod.a23b1 := ko;
kod.a34b1 := ko;
kod.a45b1 := ko;
kod.a56b1 := ko;        
kod.a67b1 := ko;        
kod.a78b1 := ko;      


dQx.b1_op=-0.035;
dQy.b1_op=-0.025;
dQpx.b1_op=15;
dQpy.b1_op=15;

twiss, file=twiss_b1;


PTC_CREATE_UNIVERSE;
PTC_CREATE_LAYOUT, MODEL=2, METHOD=6, NST=3;
ptc_twiss, NORMAL=TRUE, TRACKRDTS=TRUE, NO=4 ;
write, table=TWISSRDT , file=rdts;



stop;