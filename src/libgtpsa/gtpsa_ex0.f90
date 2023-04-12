! o-----------------------------------------------------------------------------o
! |
! | TPSA Example 0
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

program gtpsa_ex0
  use gtpsa
  implicit none

  real(c_num_t)  :: pi_6 = 3.14159265358979323846264338327950288d0/6
  type(c_ptr)    :: d, t1, t2
  real(c_num_t)  :: vec(1)

  ! descriptor for TPSA with 1 variable of order 1 (smallest construction)
  d = mad_desc_newv(1, 1_1)

  ! two TPSAs, t1 has maximum order, t2 is same as t1
  t1 = mad_tpsa_newd(d , mad_tpsa_default)
  t2 = mad_tpsa_new (t1, mad_tpsa_same)

  ! set order 0 to pi/6 and order 1 to 0 (scalar-like)
  vec = [pi_6]
  call mad_tpsa_setv (t1, 0, 1, vec)
  call mad_tpsa_print(t1, "ini"//c_eos, 0d0, 0,c_null)

  ! t2=sin(t1)
  call mad_tpsa_sin  (t1, t2)
  call mad_tpsa_print(t2, "sin"//c_eos, 1d-14, 0,c_null)
  call mad_tpsa_del  (t1); t1=c_null

  ! tpsa functions and operators support aliasing (i.e. src == dst)
  call mad_tpsa_asin (t2, t2);           ! asin(x) = -i*ln(i*x + sqrt(1-x^2))
  call mad_tpsa_print(t2, "asin"//c_eos, 1d-14, 0,c_null)
  call mad_tpsa_del(t2); t2=c_null       ! see the accuracy of asin(sin)

  ! destroy all created descriptors (optional cleanup)
  call mad_desc_cleanup()

end program gtpsa_ex0
