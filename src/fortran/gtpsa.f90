!o-----------------------------------------------------------------------------o
!|
!| GTPSA Fortran to C interface
!|
!| Methodical Accelerator Design - Copyright CERN 2016+
!| Support: http://cern.ch/mad  - mad at cern.ch
!| Authors: L. Deniau, laurent.deniau at cern.ch
!| Contrib: -
!|
!o-----------------------------------------------------------------------------o
!| You can redistribute this file and/or modify it under the terms of the GNU
!| General Public License GPLv3 (or later), as published by the Free Software
!| Foundation. This file is distributed in the hope that it will be useful, but
!| WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
!o-----------------------------------------------------------------------------o
!|
!| Notes about binding MAD C types (e.g. defs.h) to Fortran ISO C types:
!|   tpsa_t, ctpsa_t, desc_t -> type(c_ptr)          ! void* (opaque pointer)
!|   log_t                   -> logical(c_bool)
!|   idx_t, ssz_t            -> integer(c_int32_t)
!|   ord_t                   -> integer(c_signed_char)
!|   int                     -> integer(c_int)
!|   num_t                   -> real(c_double)
!|   str_t                   -> character(c_char)    ! string (*)
!|   xxx_t*                  -> ftype(c_xxx)(*)      ! pointer to array of
!|
!| Compilation:
!|   gfortran -W -Wall -Wextra -pedantic -std=f2003 -c gtpsa.f90
!o-----------------------------------------------------------------------------o

module GTPSA
  use, intrinsic :: iso_c_binding
  implicit none

  ! -- Constants ---------------------------------------------------------------

  integer(c_signed_char), parameter :: mad_tpsa_default = -1
  integer(c_signed_char), parameter :: mad_tpsa_same    = -2

  ! -- Monomials ---------------------------------------------------------------

  interface

    integer(c_int32_t) function mad_mono_str(n,mono_a,s) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n  ! monomial and string length
      integer(c_signed_char) :: mono_a(*)         ! monomial
      character(c_char), intent(in) :: s(*)       ! monomial as string "[0-9]*"
    end function mad_mono_str

    subroutine mad_mono_fill(n,mono_a,v) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n     ! monomial length
      integer(c_signed_char) :: mono_a(*)            ! monomial
      integer(c_signed_char), intent(in), value :: v ! value
    end subroutine mad_mono_fill

    subroutine mad_mono_copy(n,mono_a,mono_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*) ! src
      integer(c_signed_char) :: mono_r(*) !           ! dst
    end subroutine mad_mono_copy

    integer(c_signed_char) function mad_mono_min(n,mono_a) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomial length
      integer(c_signed_char), intent(in) :: mono_a(*) ! monomial
    end function mad_mono_min

    integer(c_signed_char) function mad_mono_max(n,mono_a) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomial length
      integer(c_signed_char), intent(in) :: mono_a(*) ! monomial
    end function mad_mono_max

    integer(c_int) function mad_mono_ord(n,mono_a) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomial length
      integer(c_signed_char), intent(in) :: mono_a(*) ! monomial
    end function mad_mono_ord

    integer(c_int) function mad_mono_eq(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_eq

    integer(c_int) function mad_mono_lt(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_lt

    integer(c_int) function mad_mono_gt(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_gt

    integer(c_int) function mad_mono_le(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_le

    integer(c_int) function mad_mono_ge(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_ge

    integer(c_int) function mad_mono_rcmp(n,mono_a,mono_b) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_rcmp

    subroutine mad_mono_add(n,mono_a,mono_b,mono_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*),  mono_b(*) ! src
      integer(c_signed_char) :: mono_r(*) !           ! dst
    end subroutine mad_mono_add

    subroutine mad_mono_sub(n,mono_a,mono_b,mono_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*),  mono_b(*) ! src
      integer(c_signed_char) :: mono_r(*) !           ! dst
    end subroutine mad_mono_sub

    subroutine mad_mono_concat(na,mono_a,nb,mono_b,mono_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb ! monomials lengths
      integer(c_signed_char), intent(in) :: mono_a(*),  mono_b(*) ! src
      integer(c_signed_char) :: mono_r(*) !           ! dst[na+nb]
    end subroutine mad_mono_concat

    subroutine mad_mono_sort(n,mono_a,idxs) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomial length
      integer(c_signed_char), intent(in) :: mono_a(*) ! src
      integer(c_int32_t) :: idxs(*) ! index lookup: a[idxs[i]] is sorted by order
    end subroutine mad_mono_sort

    subroutine mad_mono_print(n,mono_a) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! monomial length
      integer(c_signed_char), intent(in) :: mono_a(*) ! src -> stdout
    end subroutine mad_mono_print

  end interface

  ! -- Descriptor --------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a descriptor)

    type(c_ptr) function mad_desc_newn(nmv,nmo) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv         ! #map_vars
      integer(c_signed_char), intent(in), value :: nmo ! order of map_vars
    end function mad_desc_newn

    type(c_ptr) function mad_desc_newk(nmv,nmo,nk,ko,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nk     ! #map_vars, #knobs
      integer(c_signed_char), intent(in), value :: nmo, ko, dk ! order of map_vars, knobs and across knobs
    end function mad_desc_newk           ! ko <= dk <= nk*ko, dk=0 => dk=nk*ko

    type(c_ptr) function mad_desc_newm(nmv,mvar_ords) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv           ! #map_vars
      integer(c_signed_char), intent(in) :: mvar_ords(*) ! orders of map_vars
    end function mad_desc_newm

    type(c_ptr) function mad_desc_newv(nmv,mvar_ords,nv,var_ords,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nv       ! #map_vars, #vars
      integer(c_signed_char), intent(in) :: mvar_ords(*), var_ords(*) ! orders of map_vars, vars
      integer(c_signed_char), intent(in), value :: dk    ! max order across knobs
    end function mad_desc_newv  ! max(vo[nmv+1:nv]) <= dk <= ord(vo[nmv+1:nv]), dk=0 => dk=ord

    type(c_ptr) function mad_desc_newkv(nmv,mvar_ords,nkv,kvar_ords,nv,var_ords,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nkv, nv  ! #map_vars, #knobs, #vars
      integer(c_signed_char), intent(in) :: mvar_ords(*), kvar_ords(*), var_ords(*) ! orders of map_vars, knobs, vars
      integer(c_signed_char), intent(in), value :: dk    ! max order across knobs
    end function mad_desc_newkv ! max(vo[nmv+1:nv]) <= dk <= ord(vo[nmv+1:nv]), dk=0 => dk=ord

    ! -- Destructor -------------------

    subroutine mad_desc_del(desc) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
    end subroutine mad_desc_del

    subroutine mad_desc_cleanup() bind(C)
      ! global cleanup (warning: no GTSPA must still be in use!)
    end subroutine mad_desc_cleanup

    ! -- Introspection ----------------

    integer(c_signed_char) function mad_desc_maxord(desc) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
    end function mad_desc_maxord

    integer(c_int32_t) function mad_desc_maxlen(desc) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
    end function mad_desc_maxlen

    integer(c_int32_t) function mad_desc_ordlen(desc,mo) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
      integer(c_signed_char), intent(in), value :: mo ! ordlen(maxord) == maxlen
    end function mad_desc_ordlen

    integer(c_signed_char) function mad_desc_gtrunc(desc,to) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
      integer(c_signed_char), intent(in), value :: to ! global truncation order
    end function mad_desc_gtrunc

  end interface

  ! -- TPSA --------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a TPSA)

    type(c_ptr) function mad_tpsa_newd(desc,mo) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
      integer(c_signed_char), intent(in), value :: mo ! if mo > d_mo, mo = d_mo
    end function mad_tpsa_newd

    type(c_ptr) function mad_tpsa_new(tpsa,mo) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_signed_char), intent(in), value :: mo ! if mo > d_mo, mo = d_mo
    end function mad_tpsa_new

    ! -- Destructor -------------------

    subroutine mad_tpsa_del(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
    end subroutine mad_tpsa_del

    ! -- Introspection ----------------

    type(c_ptr) function mad_tpsa_desc(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
    end function mad_tpsa_desc

    integer(c_int32_t) function mad_tpsa_len(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
    end function mad_tpsa_len

    integer(c_signed_char) function mad_tpsa_ord(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
    end function mad_tpsa_ord

    ! mad_tpsa_ordv not supported by fortran

    integer(c_signed_char) function mad_tpsa_ordn(n,tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n  ! #tpsa
      type(c_ptr), intent(in) :: tpsa(*)
    end function mad_tpsa_ordn

    logical(c_bool) function mad_tpsa_is_valid(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa ! sanity check on TPSA integrity
    end function mad_tpsa_is_valid

    ! -- Initialization ---------------

    subroutine mad_tpsa_copy(tpsa,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa    ! src
      type(c_ptr), value :: tpsa_r              ! dst
    end subroutine mad_tpsa_copy

    subroutine mad_tpsa_convert(tpsa,tpsa_r,n,t2r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa     ! src
      type(c_ptr), value :: tpsa_r               ! dst
      integer(c_int32_t), intent(in), value :: n ! vector length
      integer(c_int32_t), intent(in) :: t2r(*)   ! vector of index lookup
    end subroutine mad_tpsa_convert

    subroutine mad_tpsa_clear(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
    end subroutine mad_tpsa_clear

    subroutine mad_tpsa_scalar(tpsa,v,iv,scl) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      real(c_double), intent(in), value :: v, scl  ! 0 and 1st order values
      integer(c_int32_t), intent(in), value :: iv  ! variable index (1st order)
    end subroutine mad_tpsa_scalar

    ! -- Indexing / monomials ---------

    integer(c_signed_char) function mad_tpsa_mono(tpsa,n,m,i) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n, i ! monomial length, slot index
      integer(c_signed_char) :: m(*)                ! monomial to fill
    end function mad_tpsa_mono

    integer(c_int32_t) function mad_tpsa_idxs(tpsa,n,s) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
    end function mad_tpsa_idxs

    integer(c_int32_t) function mad_tpsa_idxm(tpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end function mad_tpsa_idxm

    integer(c_int32_t) function mad_tpsa_idxsm(tpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_int32_t), intent(in) :: m(*)      ! sparse monomial (idx,ord)
    end function mad_tpsa_idxsm

    ! -- Getters ----------------------

    real(c_double) function mad_tpsa_get0(tpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
    end function mad_tpsa_get0

    real(c_double) function mad_tpsa_geti(tpsa,i) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: i ! slot index
    end function mad_tpsa_geti

    real(c_double) function mad_tpsa_gets(tpsa,n,s) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
    end function mad_tpsa_gets

    real(c_double) function mad_tpsa_getm(tpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end function mad_tpsa_getm

    real(c_double) function mad_tpsa_getsm(tpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_int32_t), intent(in) :: m(*)      ! sparse monomial (idx,ord)
    end function mad_tpsa_getsm

    subroutine mad_tpsa_getv(tpsa,i,n,v) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa
      integer(c_int32_t), intent(in), value :: i, n  ! slot index, vector length
      real(c_double) :: v(*)                         ! vector to fill
    end subroutine mad_tpsa_getv

    ! -- Setters ----------------------

    subroutine mad_tpsa_set0(tpsa,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      real(c_double), intent(in), value :: a, b ! t[0] = a*t[0]+b
    end subroutine mad_tpsa_set0

    subroutine mad_tpsa_seti(tpsa,i,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: i ! slot index
      real(c_double), intent(in), value :: a, b  ! t[i] = a*t[i]+b
    end subroutine mad_tpsa_seti

    subroutine mad_tpsa_sets(tpsa,n,s,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
      real(c_double), intent(in), value :: a, b  ! t[s] = a*t[s]+b
    end subroutine mad_tpsa_sets

    subroutine mad_tpsa_setm(tpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_signed_char), intent(in) :: m(*) ! monomial
      real(c_double), intent(in), value :: a, b  ! t[m] = a*t[m]+b
    end subroutine mad_tpsa_setm

    subroutine mad_tpsa_setsm(tpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_int32_t), intent(in) :: m(*)     ! sparse monomial (idx,ord)
      real(c_double), intent(in), value :: a, b  ! t[m] = a*t[m]+b
    end subroutine mad_tpsa_setsm

    subroutine mad_tpsa_setv(tpsa,i,n,v) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: i, n ! slot index, vector length
      real(c_double), intent(in) :: v(*)            ! vector to copy
    end subroutine mad_tpsa_setv

    ! -- Operators --------------------

    logical(c_bool) function mad_tpsa_equ(tpsa_a,tpsa_b,eps) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b
      real(c_double), intent(in), value :: eps  ! tolerance during comparison
    end function mad_tpsa_equ

    subroutine mad_tpsa_add(tpsa_a,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_add

    subroutine mad_tpsa_sub(tpsa_a,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_sub

    subroutine mad_tpsa_mul(tpsa_a,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_mul

    subroutine mad_tpsa_div(tpsa_a,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_div

    subroutine mad_tpsa_pow(tpsa_a,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_pow

    subroutine mad_tpsa_powi(tpsa_a,n,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a   ! src
      integer(c_int32_t), intent(in), value :: n ! power (integer)
      type(c_ptr), value :: tpsa_r               ! dst
    end subroutine mad_tpsa_powi

    subroutine mad_tpsa_pown(tpsa_a,v,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a   ! src
      real(c_double), intent(in), value :: v     ! power (real)
      type(c_ptr), value :: tpsa_r               ! dst
    end subroutine mad_tpsa_pown

    ! -- Functions --------------------

    subroutine mad_tpsa_abs(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a  ! src
      type(c_ptr), value :: tpsa_r              ! dst
    end subroutine mad_tpsa_abs

    real(c_double) function mad_tpsa_nrm1(tpsa_a,tpsa_b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! sum_i |a[i]-b[i]|
    end function mad_tpsa_nrm1

    real(c_double) function mad_tpsa_nrm2(tpsa_a,tpsa_b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b  ! sqrt(sum_i (a[i]-b[i])^2)
    end function mad_tpsa_nrm2

    subroutine mad_tpsa_deriv(tpsa_a,tpsa_r,iv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a     ! src
      type(c_ptr), value :: tpsa_r                 ! dst
      integer(c_int32_t), intent(in), value :: iv  ! variable index (1st order)
    end subroutine mad_tpsa_deriv

    subroutine mad_tpsa_derivm(tpsa_a,tpsa_r,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end subroutine mad_tpsa_derivm

    subroutine mad_tpsa_poisson(tpsa_a,tpsa_b,tpsa_r,nv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, tpsa_b ! src
      type(c_ptr), value :: tpsa_r                     ! dst
      integer(c_int32_t), intent(in), value :: nv    ! #variables (desc%nv if 0)
    end subroutine mad_tpsa_poisson

    subroutine mad_tpsa_taylor(tpsa_a,n,coef,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      integer(c_int32_t), intent(in), value :: n  ! vector length
      real(c_double), intent(in) :: coef(*)       ! vector of taylor coefs
    end subroutine mad_tpsa_taylor

    subroutine mad_tpsa_acc(tpsa_a,v,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! src and dst
      real(c_double), intent(in), value :: v      ! r = r + v*a (r not reset!)
    end subroutine mad_tpsa_acc

    subroutine mad_tpsa_scl(tpsa_a,v,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_double), intent(in), value :: v      ! r = v*a
    end subroutine mad_tpsa_scl

    subroutine mad_tpsa_inv(tpsa_a,v,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_double), intent(in), value :: v      ! r = v/a
    end subroutine mad_tpsa_inv

    subroutine mad_tpsa_invsqrt(tpsa_a,v,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_double), intent(in), value :: v      ! r = v/sqrt(a)
    end subroutine mad_tpsa_invsqrt

    subroutine mad_tpsa_sqrt(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sqrt

    subroutine mad_tpsa_exp(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_exp

    subroutine mad_tpsa_log(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_log

    subroutine mad_tpsa_sincos(tpsa_a,tpsa_s,tpsa_c) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_s, tpsa_c        ! dst_sin, dst_cos
    end subroutine mad_tpsa_sincos

    subroutine mad_tpsa_sin(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sin

    subroutine mad_tpsa_cos(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cos

    subroutine mad_tpsa_tan(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_tan

    subroutine mad_tpsa_cot(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cot

    subroutine mad_tpsa_sinc(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinc

    subroutine mad_tpsa_sincosh(tpsa_a,tpsa_s,tpsa_c) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_s, tpsa_c        ! dst_sin, dst_cos
    end subroutine mad_tpsa_sincosh

    subroutine mad_tpsa_sinh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinh

    subroutine mad_tpsa_cosh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cosh

    subroutine mad_tpsa_tanh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_tanh

    subroutine mad_tpsa_coth(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_coth

    subroutine mad_tpsa_sinhc(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinhc

    subroutine mad_tpsa_asin(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_asin

    subroutine mad_tpsa_acos(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acos

    subroutine mad_tpsa_atan(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_atan

    subroutine mad_tpsa_acot(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acot

    subroutine mad_tpsa_asinh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_asinh

    subroutine mad_tpsa_acosh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acosh

    subroutine mad_tpsa_atanh(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_atanh

    subroutine mad_tpsa_acoth(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acoth

    subroutine mad_tpsa_erf(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_erf

    subroutine mad_tpsa_erfc(tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_erfc

    ! -- High level functions ---------

    subroutine mad_tpsa_axpb(a,tpsa_x,b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b   ! coefs
      type(c_ptr), intent(in), value :: tpsa_x    ! src
      type(c_ptr), value :: tpsa_r                ! dst=a*x+b
    end subroutine mad_tpsa_axpb

    subroutine mad_tpsa_axpbypc(a,tpsa_x,b,tpsa_y,c,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c      ! coefs
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y  ! src
      type(c_ptr), value :: tpsa_r                      ! dst=a*x+b*y+c
    end subroutine mad_tpsa_axpbypc

    subroutine mad_tpsa_axypb(a,tpsa_x,tpsa_y,b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b         ! coefs
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y  ! src
      type(c_ptr), value :: tpsa_r                      ! dst=a*x*y+b
    end subroutine mad_tpsa_axypb

    subroutine mad_tpsa_axypbzpc(a,tpsa_x,tpsa_y,b,tpsa_z,c,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c             ! coefs
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y, tpsa_z ! src
      type(c_ptr), value :: tpsa_r                             ! dst=a*x*y+b*z+c
    end subroutine mad_tpsa_axypbzpc

    subroutine mad_tpsa_axypbvwpc(a,tpsa_x,tpsa_y,b,tpsa_u,tpsa_v,c,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c           ! coefs
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y, tpsa_u, tpsa_v ! src
      type(c_ptr), value :: tpsa_r                           ! dst=a*x*y+b*u*v+c
    end subroutine mad_tpsa_axypbvwpc

    subroutine mad_tpsa_ax2pby2pcz2(a,tpsa_x,b,tpsa_y,c,tpsa_z,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c     ! coefs
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y, tpsa_z ! src
      type(c_ptr), value :: tpsa_r                     ! dst=a*x^2+b*y^2+c*z^2
    end subroutine mad_tpsa_ax2pby2pcz2

    subroutine mad_tpsa_axpsqrtbpcx2(tpsa_x,a,b,c,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c  ! coefs
      type(c_ptr), intent(in), value :: tpsa_x      ! src
      type(c_ptr), value :: tpsa_r                  ! dst=a*x+sqrt(b+c*x^2)
    end subroutine mad_tpsa_axpsqrtbpcx2

    subroutine mad_tpsa_logaxpsqrtbpcx2(tpsa_x,a,b,c,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: a, b, c  ! coefs
      type(c_ptr), intent(in), value :: tpsa_x      ! src
      type(c_ptr), value :: tpsa_r                  ! dst=log(a*x+sqrt(b+c*x^2))
    end subroutine mad_tpsa_logaxpsqrtbpcx2

    subroutine mad_tpsa_logxdy(tpsa_x,tpsa_y,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_x, tpsa_y ! src
      type(c_ptr), value :: tpsa_r                     ! dst=log(x/y)
    end subroutine mad_tpsa_logxdy

    ! -- Maps-like functions ----------

    subroutine mad_tpsa_minv(n,tpsa_a,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst
    end subroutine mad_tpsa_minv

    subroutine mad_tpsa_pminv(n,tpsa_a,tpsa_r,select) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n      ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst
      integer(c_int32_t), intent(in) :: select(*)     ! slots to selected
    end subroutine mad_tpsa_pminv

    subroutine mad_tpsa_compose(na,tpsa_a,nb,tpsa_b,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*), tpsa_b(*) ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst[na]
    end subroutine mad_tpsa_compose

    subroutine mad_tpsa_translate(na,tpsa_a,nb,vb,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      real(c_double), intent(in) :: vb(*)             ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst[na]
    end subroutine mad_tpsa_translate

    subroutine mad_tpsa_eval(na,tpsa_a,nb,vb,vr) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      real(c_double), intent(in) :: vb(*)             ! src
      real(c_double) :: vr(*)                         ! dst[nb]
    end subroutine mad_tpsa_eval

  end interface

  ! -- CTPSA -------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a CTPSA)

    type(c_ptr) function mad_ctpsa_newd(desc,mo) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: desc
      integer(c_signed_char), intent(in), value :: mo ! if mo > d_mo, mo = d_mo
    end function mad_ctpsa_newd

    type(c_ptr) function mad_ctpsa_new(ctpsa,mo) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_signed_char), intent(in), value :: mo ! if mo > d_mo, mo = d_mo
    end function mad_ctpsa_new

    ! -- Destructor -------------------

    subroutine mad_ctpsa_del(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
    end subroutine mad_ctpsa_del

    ! -- Introspection ----------------

    type(c_ptr) function mad_ctpsa_desc(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
    end function mad_ctpsa_desc

    integer(c_int32_t) function mad_ctpsa_len(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
    end function mad_ctpsa_len

    integer(c_signed_char) function mad_ctpsa_ord(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
    end function mad_ctpsa_ord

    ! mad_ctpsa_ordv not supported by fortran

    integer(c_signed_char) function mad_ctpsa_ordn(n,ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n ! #ctpsa
      type(c_ptr), intent(in) :: ctpsa(*)
    end function mad_ctpsa_ordn

    logical(c_bool) function mad_ctpsa_is_valid(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa ! sanity check on TPSA integrity
    end function mad_ctpsa_is_valid

    ! -- Initialization ---------------

    subroutine mad_ctpsa_copy(ctpsa,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_copy

    subroutine mad_ctpsa_convert(ctpsa,ctpsa_r,n,t2r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
      integer(c_int32_t), intent(in), value :: n ! vector length
      integer(c_int32_t), intent(in) :: t2r(*)   ! vector of index lookup
    end subroutine mad_ctpsa_convert

    subroutine mad_ctpsa_clear(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
    end subroutine mad_ctpsa_clear

    subroutine mad_ctpsa_scalar(ctpsa,v,iv,scl) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      complex(c_double_complex), intent(in), value :: v, scl  ! 0 and 1st order values
      integer(c_int32_t), intent(in), value :: iv  ! variable index (1st order)
    end subroutine mad_ctpsa_scalar

    ! -- Conversion -------------------

    subroutine mad_ctpsa_real(ctpsa,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa   ! src
      type(c_ptr), value :: tpsa_r              ! dst=real(src)
    end subroutine mad_ctpsa_real

    subroutine mad_ctpsa_imag(ctpsa,tpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa   ! src
      type(c_ptr), value :: tpsa_r              ! dst=imag(src)
    end subroutine mad_ctpsa_imag

    subroutine mad_ctpsa_complex(tpsa_re,tpsa_im,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_re, tpsa_im ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=re+i*im
    end subroutine mad_ctpsa_complex

    ! -- Indexing / monomials ---------

    integer(c_signed_char) function mad_ctpsa_mono(ctpsa,n,m,i) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n, i ! monomial length, slot index
      integer(c_signed_char) :: m(*)                ! monomial to fill
    end function mad_ctpsa_mono

    integer(c_int32_t) function mad_ctpsa_idxs(ctpsa,n,s) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
    end function mad_ctpsa_idxs

    integer(c_int32_t) function mad_ctpsa_idxm(ctpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end function mad_ctpsa_idxm

    integer(c_int32_t) function mad_ctpsa_idxsm(ctpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_int32_t), intent(in) :: m(*)      ! sparse monomial (idx,ord)
    end function mad_ctpsa_idxsm

    ! -- Getters ----------------------

    complex(c_double_complex) function mad_ctpsa_get0(ctpsa) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
    end function mad_ctpsa_get0

    complex(c_double_complex) function mad_ctpsa_geti(ctpsa,i) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: i ! slot index
    end function mad_ctpsa_geti

    complex(c_double_complex) function mad_ctpsa_gets(ctpsa,n,s) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
    end function mad_ctpsa_gets

    complex(c_double_complex) function mad_ctpsa_getm(ctpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end function mad_ctpsa_getm

    complex(c_double_complex) function mad_ctpsa_getsm(ctpsa,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_int32_t), intent(in) :: m(*)      ! sparse monomial (idx,ord)
    end function mad_ctpsa_getsm

    subroutine mad_ctpsa_getv(ctpsa,i,n,v) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa
      integer(c_int32_t), intent(in), value :: i, n  ! slot index, vector length
      complex(c_double_complex) :: v(*)              ! vector to fill
    end subroutine mad_ctpsa_getv

    ! -- Setters ----------------------

    subroutine mad_ctpsa_set0(ctpsa,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      complex(c_double_complex), intent(in), value :: a, b ! ct[0] = a*ct[0]+b
    end subroutine mad_ctpsa_set0

    subroutine mad_ctpsa_seti(ctpsa,i,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      integer(c_int32_t), intent(in), value :: i            ! slot index
      complex(c_double_complex), intent(in), value :: a, b  ! ct[i] = a*ct[i]+b
    end subroutine mad_ctpsa_seti

    subroutine mad_ctpsa_sets(ctpsa,n,s,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
      complex(c_double_complex), intent(in), value :: a, b  ! ct[s] = a*ct[s]+b
    end subroutine mad_ctpsa_sets

    subroutine mad_ctpsa_setm(ctpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_signed_char), intent(in) :: m(*) ! monomial
      complex(c_double_complex), intent(in), value :: a, b ! ct[m] = a*ct[m]+b
    end subroutine mad_ctpsa_setm

    subroutine mad_ctpsa_setsm(ctpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_int32_t), intent(in) :: m(*)     ! sparse monomial (idx,ord)
      complex(c_double_complex), intent(in), value :: a, b ! ct[m] = a*ct[m]+b
    end subroutine mad_ctpsa_setsm

    subroutine mad_ctpsa_setv(ctpsa,i,n,v) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: ctpsa
      integer(c_int32_t), intent(in), value :: i, n ! slot index, vector length
      complex(c_double_complex), intent(in) :: v(*) ! vector to copy
    end subroutine mad_ctpsa_setv

    ! -- Operators --------------------

    logical(c_bool) function mad_ctpsa_equ(ctpsa_a,ctpsa_b,eps) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b
      real(c_double), intent(in), value :: eps  ! tolerance during comparison
    end function mad_ctpsa_equ

    subroutine mad_ctpsa_add(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_add

    subroutine mad_ctpsa_sub(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_sub

    subroutine mad_ctpsa_mul(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_mul

    subroutine mad_ctpsa_div(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_div

    subroutine mad_ctpsa_pow(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_pow

    subroutine mad_ctpsa_powi(ctpsa_a,n,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a   ! src
      integer(c_int32_t), intent(in), value :: n  ! power (integer)
      type(c_ptr), value :: ctpsa_r               ! dst
    end subroutine mad_ctpsa_powi

    subroutine mad_ctpsa_pown(ctpsa_a,v,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a         ! src
      complex(c_double_complex), intent(in), value :: v ! power (real)
      type(c_ptr), value :: ctpsa_r                     ! dst
    end subroutine mad_ctpsa_pown

    ! -- Operators (complex @ real, real @ complex)

    logical(c_bool) function mad_ctpsa_equt(ctpsa_a,tpsa_b,eps) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b
      real(c_double), intent(in), value :: eps  ! tolerance during comparison
    end function mad_ctpsa_equt

    subroutine mad_ctpsa_addt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_addt

    subroutine mad_ctpsa_subt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_subt

    subroutine mad_ctpsa_tsub(tpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_tsub

    subroutine mad_ctpsa_mult(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_mult

    subroutine mad_ctpsa_divt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_divt

    subroutine mad_ctpsa_tdiv(tpsa_a,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_tdiv

    subroutine mad_ctpsa_poisst(ctpsa_a,tpsa_b,ctpsa_r,nv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, tpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int32_t), intent(in), value :: nv    ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_poisst

    subroutine mad_ctpsa_tpoiss(tpsa_a,ctpsa_b,ctpsa_r,nv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: tpsa_a, ctpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int32_t), intent(in), value :: nv    ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_tpoiss

    ! -- Functions --------------------

    subroutine mad_ctpsa_abs(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_abs

    subroutine mad_ctpsa_arg(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_arg

    subroutine mad_ctpsa_conj(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_conj

    complex(c_double_complex) function mad_ctpsa_nrm1(ctpsa_a,ctpsa_b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b ! sum_i|a[i]-b[i]|
    end function mad_ctpsa_nrm1

    complex(c_double_complex) function mad_ctpsa_nrm2(ctpsa_a,ctpsa_b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a, ctpsa_b ! sqrt(sum_i(a[i]-b[i])^2)
    end function mad_ctpsa_nrm2

    subroutine mad_ctpsa_deriv(ctpsa_a,ctpsa_r,iv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      integer(c_int32_t), intent(in), value :: iv  ! variable index (1st order)
    end subroutine mad_ctpsa_deriv

    subroutine mad_ctpsa_derivm(ctpsa_a,ctpsa_r,n,m) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a   ! src
      type(c_ptr), value :: ctpsa_r               ! dst
      integer(c_int32_t), intent(in), value :: n  ! monomial length
      integer(c_signed_char), intent(in) :: m(*)  ! monomial
    end subroutine mad_ctpsa_derivm

    subroutine mad_ctpsa_poisson(ctpsa_a,ctpsa_b,ctpsa_r,nv) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a,ctpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int32_t), intent(in), value :: nv    ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_poisson

    subroutine mad_ctpsa_taylor(ctpsa_a,n,coef,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a        ! src
      type(c_ptr), value :: ctpsa_r                    ! dst
      integer(c_int32_t), intent(in), value :: n       ! vector length
      complex(c_double_complex), intent(in) :: coef(*) ! vector of taylor coefs
    end subroutine mad_ctpsa_taylor

    subroutine mad_ctpsa_acc(ctpsa_a,v,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a         ! src
      type(c_ptr), value :: ctpsa_r                     ! src and dst
      complex(c_double_complex), intent(in), value :: v ! r = r+v*a
    end subroutine mad_ctpsa_acc

    subroutine mad_ctpsa_scl(ctpsa_a,v,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a          ! src
      type(c_ptr), value :: ctpsa_r                      ! dst
      complex(c_double_complex), intent(in), value :: v  ! r = v*a
    end subroutine mad_ctpsa_scl

    subroutine mad_ctpsa_inv(ctpsa_a,v,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a          ! src
      type(c_ptr), value :: ctpsa_r                      ! dst
      complex(c_double_complex), intent(in), value :: v  ! r = v/a
    end subroutine mad_ctpsa_inv

    subroutine mad_ctpsa_invsqrt(ctpsa_a,v,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a          ! src
      type(c_ptr), value :: ctpsa_r                      ! dst
      complex(c_double_complex), intent(in), value :: v  ! r = v/sqrt(a)
    end subroutine mad_ctpsa_invsqrt

    subroutine mad_ctpsa_sqrt(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a     ! src
      type(c_ptr), value :: ctpsa_r                 ! dst
    end subroutine mad_ctpsa_sqrt

    subroutine mad_ctpsa_exp(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_exp

    subroutine mad_ctpsa_log(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_log

    subroutine mad_ctpsa_sincos(ctpsa_a,ctpsa_s,ctpsa_c) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_s, ctpsa_c       ! dst_sin, dst_cos
    end subroutine mad_ctpsa_sincos

    subroutine mad_ctpsa_sin(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sin

    subroutine mad_ctpsa_cos(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cos

    subroutine mad_ctpsa_tan(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_tan

    subroutine mad_ctpsa_cot(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cot

    subroutine mad_ctpsa_sinc(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinc

    subroutine mad_ctpsa_sincosh(ctpsa_a,ctpsa_s,ctpsa_c) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_s, ctpsa_c       ! dst_sin, dst_cos
    end subroutine mad_ctpsa_sincosh

    subroutine mad_ctpsa_sinh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinh

    subroutine mad_ctpsa_cosh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cosh

    subroutine mad_ctpsa_tanh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_tanh

    subroutine mad_ctpsa_coth(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_coth

    subroutine mad_ctpsa_sinhc(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinhc

    subroutine mad_ctpsa_asin(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_asin

    subroutine mad_ctpsa_acos(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acos

    subroutine mad_ctpsa_atan(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_atan

    subroutine mad_ctpsa_acot(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acot

    subroutine mad_ctpsa_asinh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_asinh

    subroutine mad_ctpsa_acosh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acosh

    subroutine mad_ctpsa_atanh(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_atanh

    subroutine mad_ctpsa_acoth(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acoth

    subroutine mad_ctpsa_erf(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_erf

    subroutine mad_ctpsa_erfc(ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_erfc

    ! -- High level functions ---------

    subroutine mad_ctpsa_axpb(a,ctpsa_x,b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x            ! src
      type(c_ptr), value :: ctpsa_r                        ! dst=a*x+b
    end subroutine mad_ctpsa_axpb

    subroutine mad_ctpsa_axpbypc(a,ctpsa_x,b,ctpsa_y,c,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y      ! src
      type(c_ptr), value :: ctpsa_r                           ! dst=a*x+b*y+c
    end subroutine mad_ctpsa_axpbypc

    subroutine mad_ctpsa_axypb(a,ctpsa_x,ctpsa_y,b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b  ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y    ! src
      type(c_ptr), value :: ctpsa_r                         ! dst=a*x*y+b
    end subroutine mad_ctpsa_axypb

    subroutine mad_ctpsa_axypbzpc(a,ctpsa_x,ctpsa_y,b,ctpsa_z,c,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c     ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y, ctpsa_z ! src
      type(c_ptr), value :: ctpsa_r                            ! dst=a*x*y+b*z+c
    end subroutine mad_ctpsa_axypbzpc

    subroutine mad_ctpsa_axypbvwpc(a,ctpsa_x,ctpsa_y,b,ctpsa_u,ctpsa_v,c,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c            ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y, ctpsa_u, ctpsa_v ! src
      type(c_ptr), value :: ctpsa_r                          ! dst=a*x*y+b*u*v+c
    end subroutine mad_ctpsa_axypbvwpc

    subroutine mad_ctpsa_ax2pby2pcz2(a,ctpsa_x,b,ctpsa_y,c,ctpsa_z,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c     ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y, ctpsa_z ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=a*x^2+b*y^2+c*z^2
    end subroutine mad_ctpsa_ax2pby2pcz2

    subroutine mad_ctpsa_axpsqrtbpcx2(ctpsa_x,a,b,c,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x               ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=a*x+sqrt(b+c*x^2)
    end subroutine mad_ctpsa_axpsqrtbpcx2

    subroutine mad_ctpsa_logaxpsqrtbpcx2(ctpsa_x,a,b,c,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      complex(c_double_complex), intent(in), value :: a, b, c ! coefs
      type(c_ptr), intent(in), value :: ctpsa_x               ! src
      type(c_ptr), value :: ctpsa_r                 ! dst=log(a*x+sqrt(b+c*x^2))
    end subroutine mad_ctpsa_logaxpsqrtbpcx2

    subroutine mad_ctpsa_logxdy(ctpsa_x,ctpsa_y,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: ctpsa_x, ctpsa_y ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=log(x/y)
    end subroutine mad_ctpsa_logxdy

    ! -- Maps-like functions ----------

    subroutine mad_ctpsa_minv(n,ctpsa_a,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n       ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)            ! src
      type(c_ptr) :: ctpsa_r(*)                        ! dst
    end subroutine mad_ctpsa_minv

    subroutine mad_ctpsa_pminv(n,ctpsa_a,ctpsa_r,select) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: n       ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)            ! src
      type(c_ptr) :: ctpsa_r(*)                        ! dst
      integer(c_int32_t), intent(in) :: select(*)      ! slots to selected
    end subroutine mad_ctpsa_pminv

    subroutine mad_ctpsa_compose(na,ctpsa_a,nb,ctpsa_b,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb   ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*), ctpsa_b(*) ! src
      type(c_ptr) :: ctpsa_r(*)                         ! dst[na]
    end subroutine mad_ctpsa_compose

    subroutine mad_ctpsa_translate(na,ctpsa_a,nb,vb,ctpsa_r) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb  ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)            ! src
      complex(c_double_complex), intent(in) :: vb(*)   ! src
      type(c_ptr) :: ctpsa_r(*)                        ! dst[na]
    end subroutine mad_ctpsa_translate

    subroutine mad_ctpsa_eval(na,ctpsa_a,nb,vb,vr) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), intent(in), value :: na, nb  ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)            ! src
      complex(c_double_complex), intent(in) :: vb(*)   ! src
      complex(c_double_complex) :: vr(*)               ! dst[nb]
    end subroutine mad_ctpsa_eval

  end interface

end module GTPSA
