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

module GTPSA
  use, intrinsic :: iso_c_binding
  implicit none

  ! -- Constants ---------------------------------------------------------------

  integer(c_signed_char), parameter :: mad_tpsa_default = -1
  integer(c_signed_char), parameter :: mad_tpsa_same    = -2

  ! -- Descriptor --------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a descriptor)

    type(c_ptr) function mad_desc_newn(nmv,nmo) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv
      integer(c_signed_char), intent(in), value :: nmo
    end function mad_desc_newn

    type(c_ptr) function mad_desc_newk(nmv,nmo,nk,ko,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nk
      integer(c_signed_char), intent(in), value :: nmo, ko, dk
    end function mad_desc_newk

    type(c_ptr) function mad_desc_newm(nmv,mvar_ords) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv
      integer(c_signed_char), intent(in) :: mvar_ords(*)
    end function mad_desc_newm

    type(c_ptr) function mad_desc_newv(nmv,mvar_ords,nv,var_ords,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nv
      integer(c_signed_char), intent(in) :: mvar_ords(*), var_ords(*)
      integer(c_signed_char), intent(in), value :: dk
    end function mad_desc_newv

    type(c_ptr) function mad_desc_newkv(nmv,mvar_ords,nkv,kvar_ords,nv,var_ords,dk) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: nmv, nkv, nv
      integer(c_signed_char), intent(in) :: mvar_ords(*), kvar_ords(*), var_ords(*)
      integer(c_signed_char), intent(in), value :: dk
    end function mad_desc_newkv

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
      real(c_double), intent(in), value :: a, b ! tpsa[0] = a*tpsa[0]+b
    end subroutine mad_tpsa_set0

    subroutine mad_tpsa_seti(tpsa,i,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: i ! slot index
      real(c_double), intent(in), value :: a, b  ! tpsa[i] = a*tpsa[i]+b
    end subroutine mad_tpsa_seti

    subroutine mad_tpsa_sets(tpsa,n,s,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
      real(c_double), intent(in), value :: a, b  ! tpsa[s] = a*tpsa[s]+b
    end subroutine mad_tpsa_sets

    subroutine mad_tpsa_setm(tpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_signed_char), intent(in) :: m(*) ! monomial
      real(c_double), intent(in), value :: a, b  ! tpsa[m] = a*tpsa[m]+b
    end subroutine mad_tpsa_setm

    subroutine mad_tpsa_setsm(tpsa,n,m,a,b) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: tpsa
      integer(c_int32_t), intent(in), value :: n ! monomial length
      integer(c_int32_t), intent(in) :: m(*)     ! sparse monomial (idx,ord)
      real(c_double), intent(in), value :: a, b  ! tpsa[m] = a*tpsa[m]+b
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
      type(c_ptr), intent(in), value :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      integer(c_int32_t), intent(in), value :: nv ! #variables (desc%nv if 0)
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

  end interface

  ! -- CTPSA -------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a ctpsa)

  end interface

end module GTPSA
