! o-----------------------------------------------------------------------------o
! |
! | TPSA Example 8
! |
! | Methodical Accelerator Design - Copyright (c) 2016+
! | Support: http://cern.ch/mad  - mad at cern.ch
! | Authors: L. Deniau, laurent.deniau at cern.ch
! | Contrib: -
! |
! o-----------------------------------------------------------------------------o
! | You can redistribute this file and/or modify it under the terms of the GNU
! | General Public License GPLv3 (or later), as published by the Free Software
! | Foundation. This file is distributed in the hope that it will be useful, but
! | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
! o-----------------------------------------------------------------------------o
!
!  WARNING: this file is not part of MAD, it must be used as an example with the
!           standalone library, see mad/src/REAME.GTPSA for more info.

program gtpsa_ex8
  use gtpsa
  implicit none

  type(c_ptr) :: d, d0, t1, t2, t3

  ! descriptor for TPSA with 6 variables of order 2 (max)
  d=mad_desc_newv(6, 2)

  ! d0=mad_desc_newv(6, 0)
  ! -> error: mad_desc.c:1230: : invalid maximum order: 0 (0<?<=63)

  ! create three TPSAs of order 0 (as scalars) and set values
  t1=mad_tpsa_newd(d, 0) ; call mad_tpsa_set0(t1,0d0,1.5d0);
  t2=mad_tpsa_newd(d, 0) ; call mad_tpsa_set0(t2,0d0,2.0d0);

   call mad_tpsa_print(t1,"T1"//c_eos,0d0,0,c_null);
   call mad_tpsa_print(t2,"T2"//c_eos,0d0,0,c_null);

   ! add the scalars
   call mad_tpsa_add(t1,t2,t3);
   call mad_tpsa_print(t3,"T1+T2"//c_eos,0d0,0,c_null);

   ! multiply the scalars
   call mad_tpsa_mul(t1,t2,t3);
   call mad_tpsa_print(t3,"T1*T2"//c_eos,0d0,0,c_null);

   ! log of the scalar
   call mad_tpsa_log(t3,t3);
   call mad_tpsa_print(t3,"log(T1*T2)"//c_eos,0d0,0,c_null);

   ! destroy the three TPSAs
   call mad_tpsa_del(t1) ; t1=c_null
   call mad_tpsa_del(t2) ; t2=c_null
   call mad_tpsa_del(t3) ; t3=c_null

  ! destroy all created descriptors (optional cleanup)
  call mad_desc_cleanup();

end program gtpsa_ex8
