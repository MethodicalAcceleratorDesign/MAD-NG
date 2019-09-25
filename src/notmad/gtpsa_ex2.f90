! o-----------------------------------------------------------------------------o
! |
! | TPSA Example 2
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

program gtpsa_ex2
  use gtpsa
  implicit none

  double precision  :: pi_6 = 3.14159265358979323846264338327950288d0/6
  character(c_char) :: eos = c_null_char ! end of C string

  type(c_ptr) :: d, t1, t2
  real(c_double) :: vec(1:12)
  !integer(c_signed_char) :: no=3, ko=2, dk=2

  ! descriptor for TPSA with 6 variables of order 3 and
  !                          5 knobs of order 2 with cross-terms of order 2
  d=mad_desc_newk(6, 3_1, 5, 2_1, 2_1)

  ! two TPSAs, t2 is same as t1
  t1=mad_tpsa_newd(d , mad_tpsa_default)
  t2=mad_tpsa_new (t1, mad_tpsa_same)

  ! set order 0 and 1 (quick and dirty!)
  vec = [pi_6, 1d0,1d0,1d0,1d0,1d0,1d0, 1d0,1d0,1d0,1d0,1d0]
  call mad_tpsa_setv (t1, 0, 1+6+5, vec);
  call mad_tpsa_print(t1, "ini"//eos, 0d0, c_null_ptr);

  ! t2=sin(t1)
  call mad_tpsa_sin  (t1, t2)
  call mad_tpsa_print(t2, "sin"//eos, 0d0, c_null_ptr);
  call mad_tpsa_del  (t1)

  ! tpsa functions and operators support aliasing (i.e. src == dst)
  call mad_tpsa_asin (t2, t2);           ! asin(x) = -i*ln(i*x + sqrt(1-x^2))
  call mad_tpsa_print(t2, "asin"//eos, 0d0, c_null_ptr);
  call mad_tpsa_del  (t2);                ! see the accuracy of asin(sin)

  ! destroy all created descriptors (optional cleanup)
  call mad_desc_cleanup();

end program gtpsa_ex2
