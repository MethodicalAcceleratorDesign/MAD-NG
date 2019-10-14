!o-----------------------------------------------------------------------------o
!|
!| GTPSA Fortran to C interface
!|
!| Methodical Accelerator Design - Copyright (c) 2016+
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
!| Notes about optional arguments:
!|   By convention, identifiers of formal parameters ending with an underscore
!|   (i.e. _) have optional/default values. Hence optional integers and reals
!|   can have the value 0 (or 0_1) and 0d0, and pointers can have the value
!|   c_null (i.e. c_null_ptr).
!|
!| Notes about strings:
!|   All strings are expected to be C null terminated strings. An easy way to
!|   build null terminated strings is to concatenate c_eos (i.e. c_null_char) at
!|   the end, e.g. "example of C-style string in fortran"//c_eos.
!|
!| Notes about binding MAD C types (e.g. defs.h) to Fortran ISO C types:
!|   tpsa_t, ctpsa_t, desc_t -> type(c_ptr)          ! void* (opaque pointer)
!|   log_t                   -> logical(c_bool)
!|   idx_t, ssz_t            -> integer(c_int32_t)   ! literals: 5
!|   ord_t                   -> integer(c_ord_t)     ! literals: 5_1
!|   int                     -> integer(c_int)       ! literals: 5
!|   num_t                   -> real(c_double)       ! literals: 5d0
!|   cnum_t                  -> complex(c_double_complex)
!|   str_t                   -> character(c_char)    ! literals: "hi"//c_eos
!|   xxx_t*                  -> ftype(c_xxx)(*)      ! pointer to array of
!|   optional pointer (NULL) -> c_null_ptr           ! pointer name ending by _
!|
!| Compilation:
!|   note that f2003-2008 don't support interoperability of optional for bind(c)
!|   gfortran -W -Wall -Wextra -pedantic -c gtpsa.f90
!|   gfortran -W -Wall -Wextra -pedantic -std=f2018 -c gtpsa.f90
!o-----------------------------------------------------------------------------o

module GTPSA
  use, intrinsic :: iso_c_binding
  implicit none

  ! ----------------------------------------------------------------------------
  ! -- Types sizes -------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  integer, parameter :: c_ord_t  = c_signed_char
  integer, parameter :: c_idx_t  = c_int32_t
  integer, parameter :: c_ssz_t  = c_int32_t
  integer, parameter :: c_num_t  = c_double
  integer, parameter :: c_cnum_t = c_double_complex

  ! ----------------------------------------------------------------------------
  ! -- Constants ---------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  character(c_char), parameter :: c_eos  = c_null_char ! C-style end of string
  type(c_ptr)      , parameter :: c_null = c_null_ptr  ! C-style NULL pointer

  integer(c_ord_t) , parameter :: mad_tpsa_default = -1
  integer(c_ord_t) , parameter :: mad_tpsa_same    = -2

  ! ----------------------------------------------------------------------------
  ! -- Monomials ---------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  interface

    function mad_mono_str(n,mono_a,s) result(size) bind(C)
      import ; implicit none
      integer(c_ssz_t) :: size                 ! adjusted size n if "\0" found
      integer(c_ssz_t), value, intent(in) :: n ! monomial and string length
      integer(c_ord_t) :: mono_a(*)            ! monomial
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_mono_str

    subroutine mad_mono_fill(n,mono_a,v) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: mono_a(*)            ! monomial
      integer(c_ord_t), value, intent(in) :: v ! value
    end subroutine mad_mono_fill

    subroutine mad_mono_copy(n,mono_a,mono_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n  ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*) ! src
      integer(c_ord_t) :: mono_r(*) !           ! dst
    end subroutine mad_mono_copy

    subroutine mad_mono_rcopy(n,mono_a,mono_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n  ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*) ! src
      integer(c_ord_t) :: mono_r(*) !           ! dst
    end subroutine mad_mono_rcopy

    function mad_mono_min(n,mono_a) result(min) bind(C)
      import ; implicit none
      integer(c_ord_t) :: min                    ! min order
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: mono_a(*)  ! monomial
    end function mad_mono_min

    function mad_mono_max(n,mono_a) result(max) bind(C)
      import ; implicit none
      integer(c_ord_t) :: max                    ! max order
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: mono_a(*)  ! monomial
    end function mad_mono_max

    function mad_mono_ord(n,mono_a) result(ord) bind(C)
      import ; implicit none
      integer(c_int) :: ord                      ! order of monomial (sum)
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: mono_a(*)  ! monomial
    end function mad_mono_ord

    function mad_mono_eq(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                               ! true or false
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_eq

    function mad_mono_lt(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                               ! true or false
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_lt

    function mad_mono_gt(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                               ! true or false
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_gt

    function mad_mono_le(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                               ! true or false
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_le

    function mad_mono_ge(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                               ! true or false
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_ge

    function mad_mono_cmp(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none
      integer(c_int) :: ret                                ! first a[i]-b[i] != 0
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_cmp

    function mad_mono_rcmp(n,mono_a,mono_b) result(ret) bind(C)
      import ; implicit none                               ! compare from end
      integer(c_int) :: ret                                ! first a[i]-b[i] != 0
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! monomials
    end function mad_mono_rcmp

    subroutine mad_mono_add(n,mono_a,mono_b,mono_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! src
      integer(c_ord_t) :: mono_r(*)                        ! dst
    end subroutine mad_mono_add

    subroutine mad_mono_sub(n,mono_a,mono_b,mono_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n             ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! src
      integer(c_ord_t) :: mono_r(*)                        ! dst
    end subroutine mad_mono_sub

    subroutine mad_mono_cat(na,mono_a,nb,mono_b,mono_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb        ! monomials lengths
      integer(c_ord_t), intent(in) :: mono_a(*), mono_b(*) ! src
      integer(c_ord_t) :: mono_r(*)                        ! dst[na+nb]
    end subroutine mad_mono_cat

    subroutine mad_mono_sort(n,mono_a,idxs) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: mono_a(*)  ! src
      integer(c_idx_t) :: idxs(*)  ! index lookup: a[idxs[i]] is sorted by order
    end subroutine mad_mono_sort

    subroutine mad_mono_print(n,mono_a) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: mono_a(*)  ! src -> stdout
    end subroutine mad_mono_print

  end interface

  ! ----------------------------------------------------------------------------
  ! -- Descriptors -------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a descriptor)

    function mad_desc_newn(nv,mo_) result(desc) bind(C)
      import ; implicit none
      type(c_ptr) :: desc                          ! descriptor
      integer(c_int), value, intent(in) :: nv      ! #vars
      integer(c_ord_t), value, intent(in) :: mo_   ! order of tpsa, mo=max(1,mo_)
    end function mad_desc_newn

    function mad_desc_newk(nv,mo_,nk,ko_) result(desc) bind(C)
      ! if nk == 0, same as mad_desc_newn, otherwise
      ! mo = max(1,mo_) and ko = ko_ ? min(mo,ko_) : mo
      import ; implicit none
      type(c_ptr) :: desc                             ! descriptor
      integer(c_int), value, intent(in) :: nv, nk     ! #vars, #knobs (right part of vars)
      integer(c_ord_t), value, intent(in) :: mo_, ko_ ! order of tpsa and knobs
    end function mad_desc_newk

    function mad_desc_newv(nv,vo,nk,ko_) result(desc) bind(C)
      ! mo = max(vo[0:nv-1])
      ! ko = nk>0 ? min(mo, max(ko_, max( vo[nv-nk:nv-1] ))) : mo
      import ; implicit none
      type(c_ptr) :: desc                         ! descriptor
      integer(c_int), value, intent(in) :: nv, nk ! #vars, #knobs (i.e. mo=max(vo))
      integer(c_ord_t), value, intent(in) :: ko_  ! max order of knobs
      integer(c_ord_t), intent(in) :: vo(*)       ! orders of vars, (mvars and knobs)
    end function mad_desc_newv

    ! -- Destructor -------------------

    subroutine mad_desc_del(desc) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: desc      ! descriptor to delete
    end subroutine mad_desc_del

    subroutine mad_desc_cleanup() bind(C)
      ! global cleanup (warning: no GTSPA must still be in use!)
    end subroutine mad_desc_cleanup

    ! -- Introspection ----------------

    function mad_desc_nvmok(desc,mo_,nk_,ko_) result(nv) bind(C)
      import ; implicit none
      integer(c_int) :: nv                                ! #variables
      type(c_ptr), value, intent(in) :: desc              ! descriptor
      integer(c_int), optional, intent(out) :: nk_        ! #knobs
      integer(c_ord_t), optional, intent(out) :: mo_, ko_ ! tpsa order, knobs orders
    end function mad_desc_nvmok

    function mad_desc_getvo(desc,nv,vo_) result(mo) bind(C)
      import ; implicit none
      integer(c_ord_t) :: mo                   ! tpsa max order
      type(c_ptr), value, intent(in) :: desc   ! descriptor
      integer(c_int), value, intent(in) :: nv  ! #variables, vo_[1..nv]
      integer(c_ord_t), intent(out) :: vo_(*)  ! orders to be filled if provided
    end function mad_desc_getvo

    function mad_desc_maxord(desc) result(mo) bind(C)
      import ; implicit none
      integer(c_ord_t) :: mo                   ! tpsa max order
      type(c_ptr), value, intent(in) :: desc   ! descriptor
    end function mad_desc_maxord

    function mad_desc_maxlen(desc) result(mlen) bind(C)
      import ; implicit none
      integer(c_ssz_t) :: mlen                 ! #monomials in 0..maxorder
      type(c_ptr), value, intent(in) :: desc   ! descriptor
    end function mad_desc_maxlen

    function mad_desc_ordlen(desc,mo) result(olen) bind(C)
      import ; implicit none
      integer(c_ssz_t) :: olen                  ! #monomials in 0..order
      type(c_ptr), value, intent(in) :: desc    ! descriptor
      integer(c_ord_t), value, intent(in) :: mo ! ordlen(maxord) == maxlen
    end function mad_desc_ordlen

    function mad_desc_gtrunc(desc,to) result(oldto) bind(C)
      import ; implicit none
      integer(c_ord_t) :: oldto                 ! prev. global truncation order
      type(c_ptr), value, intent(in) :: desc    ! descriptor
      integer(c_ord_t), value, intent(in) :: to ! new global truncation order
    end function mad_desc_gtrunc

    ! -- indexes / monomials (return idx_t = -1 if invalid)

    function mad_desc_isvalids(desc,n,s) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                   ! true or false
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_desc_isvalids

    function mad_desc_isvalidm(desc,n,m) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                   ! true or false
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(in) :: m(*)     ! monomial
    end function mad_desc_isvalidm

    function mad_desc_isvalidsm(desc,n,m) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                   ! true or false
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_int), intent(in) :: m(*)       ! sparse monomial (idx,ord)
    end function mad_desc_isvalidsm

    function mad_desc_idxs(desc,n,s) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_desc_idxs

    function mad_desc_idxm(desc,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(in) :: m(*)     ! monomial
    end function mad_desc_idxm

    function mad_desc_idxsm(desc,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_int), intent(in) :: m(*)       ! sparse monomial (idx,ord)
    end function mad_desc_idxsm

    function mad_desc_nxtbyvar(desc,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(out) :: m(*)     ! monomial
    end function mad_desc_nxtbyvar

    function mad_desc_nxtbyord(desc,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(out) :: m(*)    ! monomial
    end function mad_desc_nxtbyord

    function mad_desc_mono(desc,n,m_,i) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                  ! monomial order
      type(c_ptr), value, intent(in) :: desc   !
      integer(c_idx_t), value, intent(in) :: i ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: m_(*)                ! monomial to fill (if provided)
    end function mad_desc_mono

  end interface

  ! ----------------------------------------------------------------------------
  ! -- Real GTPSA --------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a TPSA)

    function mad_tpsa_newd(desc,mo) result(newtpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: newtpsa                    ! new tpsa
      type(c_ptr), value, intent(in) :: desc    ! descriptor
      integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
    end function mad_tpsa_newd

    function mad_tpsa_new(tpsa,mo) result(newtpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: newtpsa                    ! new tpsa
      type(c_ptr), value, intent(in) :: tpsa    ! (reference) tpsa
      integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
    end function mad_tpsa_new

    ! -- Destructor -------------------

    subroutine mad_tpsa_del(tpsa) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa    ! tpsa to delete
    end subroutine mad_tpsa_del

    ! -- Introspection ----------------

    function mad_tpsa_desc(tpsa) result(desc) bind(C)
      import ; implicit none
      type(c_ptr) :: desc                       ! tpsa descriptor
      type(c_ptr), value, intent(in) :: tpsa
    end function mad_tpsa_desc

    function mad_tpsa_uid(tpsa, uid_) result(uid) bind(C)
      import ; implicit none
      integer(c_int32_t) :: uid                     ! tpsa special uid
      type(c_ptr), value, intent(in) :: tpsa
      integer(c_int32_t), value, intent(in) :: uid_ ! tpsa new uid if != 0
    end function mad_tpsa_uid

    function mad_tpsa_len(tpsa) result(len) bind(C)
      import ; implicit none
      integer(c_ssz_t) :: len                   ! #monomials in tpsa
      type(c_ptr), value, intent(in) :: tpsa
    end function mad_tpsa_len

    function mad_tpsa_ord(tpsa) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                   ! tpsa order
      type(c_ptr), value, intent(in) :: tpsa
    end function mad_tpsa_ord

    ! mad_tpsa_ordv not supported by fortran

    function mad_tpsa_ordn(n,tpsa) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                   ! max of all tpsas order
      integer(c_ssz_t), value, intent(in) :: n  ! #tpsa
      type(c_ptr), intent(in) :: tpsa(*)
    end function mad_tpsa_ordn

    function mad_tpsa_isvalid(tpsa) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                 ! true or false
      type(c_ptr), value, intent(in) :: tpsa ! sanity check on TPSA integrity
    end function mad_tpsa_isvalid

    ! -- Initialization ---------------

    subroutine mad_tpsa_copy(tpsa,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa    ! src
      type(c_ptr), value :: tpsa_r              ! dst
    end subroutine mad_tpsa_copy

    subroutine mad_tpsa_getord(tpsa,tpsa_r,ord) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa     ! src
      type(c_ptr), value :: tpsa_r               ! dst
      integer(c_ord_t), value, intent(in) :: ord ! order to retrieve
    end subroutine mad_tpsa_getord

    subroutine mad_tpsa_cutord(tpsa,tpsa_r,ord) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa    ! src
      type(c_ptr), value :: tpsa_r              ! dst
      integer(c_int), value, intent(in) :: ord  ! cut order: 0..-ord or ord..mo
    end subroutine mad_tpsa_cutord

    subroutine mad_tpsa_convert(tpsa,tpsa_r,n,t2r_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa    ! src
      type(c_ptr), value :: tpsa_r              ! dst
      integer(c_ssz_t), value, intent(in) :: n  ! vector length
      integer(c_idx_t), intent(in) :: t2r_(*)   ! vector of index lookup
    end subroutine mad_tpsa_convert

    subroutine mad_tpsa_setvar(tpsa,v,iv_,scl_) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      real(c_num_t), value, intent(in) :: v, scl_ ! 0th and 1st order values
      integer(c_idx_t), value, intent(in) :: iv_  ! variable index
    end subroutine mad_tpsa_setvar                ! equiv. to set0 if iv=0

    subroutine mad_tpsa_clear(tpsa) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa                  ! clear tpsa (reset to 0)
    end subroutine mad_tpsa_clear

    ! -- Indexing / monomials (return idx_t = -1 if invalid)

    function mad_tpsa_mono(tpsa,n,m_,i) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                  ! monomial order
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_idx_t), value, intent(in) :: i ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: m_(*)                ! monomial to fill (if provided)
    end function mad_tpsa_mono

    function mad_tpsa_idxs(tpsa,n,s) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_ssz_t), value, intent(in) :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_tpsa_idxs

    function mad_tpsa_idxm(tpsa,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(in) :: m(*)     ! monomial
    end function mad_tpsa_idxm

    function mad_tpsa_idxsm(tpsa,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index or -1
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_int), intent(in) :: m(*)       ! sparse monomial (idx,ord)
    end function mad_tpsa_idxsm

    function mad_tpsa_cycle(tpsa,n,m_,i,v_) result(idx) bind(C)
      import ; implicit none                   ! scan for non-zero coefs starting at i
      integer(c_idx_t) :: idx                  ! next index to start searching or -1
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_idx_t), value, intent(in) :: i ! index to start searching
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: m_(*)                ! monomial to fill (if provided)
      real(c_num_t), intent(out) :: v_         ! coeff to fill (if provided)
    end function mad_tpsa_cycle

    ! -- Getters ----------------------

    function mad_tpsa_get0(tpsa) result(val) bind(C)
      import ; implicit none
      real(c_num_t) :: val                     ! value at order 0 (index 0)
      type(c_ptr), value, intent(in) :: tpsa   !
    end function mad_tpsa_get0

    function mad_tpsa_geti(tpsa,i) result(val) bind(C)
      import ; implicit none
      real(c_num_t) :: val                     ! value at index i
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_idx_t), value, intent(in) :: i ! slot index (must be valid)
    end function mad_tpsa_geti

    function mad_tpsa_gets(tpsa,n,s) result(val) bind(C)
      import ; implicit none
      real(c_num_t) :: val                     ! value at string monomial
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_ssz_t), value, intent(in) :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_tpsa_gets

    function mad_tpsa_getm(tpsa,n,m) result(val) bind(C)
      import ; implicit none
      real(c_num_t) :: val                     ! value at monomial
      type(c_ptr), value, intent(in) :: tpsa   !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(in) :: m(*)     ! monomial
    end function mad_tpsa_getm

    function mad_tpsa_getsm(tpsa,n,m) result(val) bind(C)
      import ; implicit none
      real(c_num_t) :: val                      ! value at sparse monomial
      type(c_ptr), value, intent(in) :: tpsa
      integer(c_ssz_t), value, intent(in) :: n  ! monomial length
      integer(c_int), intent(in) :: m(*)        ! sparse monomial (idx,ord)
    end function mad_tpsa_getsm

    subroutine mad_tpsa_getv(tpsa,i,n,v) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa
      integer(c_idx_t), value, intent(in) :: i  ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n  ! vector length
      real(c_num_t) :: v(*)                     ! vector to fill
    end subroutine mad_tpsa_getv

    ! -- Setters ----------------------

    subroutine mad_tpsa_set0(tpsa,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      real(c_num_t), value, intent(in) :: a, b   ! t[0] = a*t[0]+b
    end subroutine mad_tpsa_set0

    subroutine mad_tpsa_seti(tpsa,i,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      integer(c_idx_t), value, intent(in) :: i   ! slot index (must be valid)
      real(c_num_t), value, intent(in) :: a, b   ! t[i] = a*t[i]+b
    end subroutine mad_tpsa_seti

    subroutine mad_tpsa_sets(tpsa,n,s,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      integer(c_ssz_t), value, intent(in) :: n   ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)      ! monomial as string "[0-9]*"
      real(c_num_t), value, intent(in) :: a, b   ! t[s] = a*t[s]+b
    end subroutine mad_tpsa_sets

    subroutine mad_tpsa_setm(tpsa,n,m,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: m(*)       ! monomial
      real(c_num_t), value, intent(in) :: a, b   ! t[m] = a*t[m]+b
    end subroutine mad_tpsa_setm

    subroutine mad_tpsa_setsm(tpsa,n,m,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_int), intent(in) :: m(*)         ! sparse monomial (idx,ord)
      real(c_num_t), value, intent(in) :: a, b   ! t[m] = a*t[m]+b
    end subroutine mad_tpsa_setsm

    subroutine mad_tpsa_setv(tpsa,i,n,v) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa
      integer(c_idx_t), value, intent(in) :: i   ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n   ! vector length
      real(c_num_t), intent(in) :: v(*)          ! vector to copy
    end subroutine mad_tpsa_setv

    ! -- Operators --------------------

    function mad_tpsa_equ(tpsa_a,tpsa_b,eps_) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                    ! true or false
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b
      real(c_num_t), value, intent(in) :: eps_  ! tolerance during comparison
    end function mad_tpsa_equ

    subroutine mad_tpsa_add(tpsa_a,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_add

    subroutine mad_tpsa_sub(tpsa_a,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_sub

    subroutine mad_tpsa_mul(tpsa_a,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_mul

    subroutine mad_tpsa_div(tpsa_a,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_div

    subroutine mad_tpsa_pow(tpsa_a,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_pow

    subroutine mad_tpsa_powi(tpsa_a,n,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a          ! src
      integer(c_int), value, intent(in) :: n            ! power (integer)
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_powi

    subroutine mad_tpsa_pown(tpsa_a,v,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a          ! src
      real(c_num_t), value, intent(in) :: v             ! power (real)
      type(c_ptr), value :: tpsa_r                      ! dst
    end subroutine mad_tpsa_pown

    ! -- Functions --------------------

    subroutine mad_tpsa_abs(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a   ! src
      type(c_ptr), value :: tpsa_r               ! dst=|src|
    end subroutine mad_tpsa_abs

    function mad_tpsa_nrm1(tpsa_a,tpsa_b_) result(nrm1) bind(C)
      import ; implicit none
      real(c_num_t) :: nrm1                      ! sum_i |a[i]-b_[i]|
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b_
    end function mad_tpsa_nrm1

    function mad_tpsa_nrm2(tpsa_a,tpsa_b_) result(nrm2) bind(C)
      import ; implicit none
      real(c_num_t) :: nrm2                      ! sqrt(sum_i (a[i]-b_[i])^2)
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b_
    end function mad_tpsa_nrm2

    subroutine mad_tpsa_deriv(tpsa_a,tpsa_r,iv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a   ! src
      type(c_ptr), value :: tpsa_r               ! dst
      integer(c_int), value, intent(in) :: iv    ! variable number (1st order)
    end subroutine mad_tpsa_deriv

    subroutine mad_tpsa_derivm(tpsa_a,tpsa_r,n,m) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a   ! src
      type(c_ptr), value :: tpsa_r               ! dst
      integer(c_ssz_t), value, intent(in) :: n   ! monomial length
      integer(c_ord_t), intent(in) :: m(*)       ! monomial
    end subroutine mad_tpsa_derivm

    subroutine mad_tpsa_poisson(tpsa_a,tpsa_b,tpsa_r,nv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, tpsa_b ! src
      type(c_ptr), value :: tpsa_r                     ! dst
      integer(c_int), value, intent(in) :: nv          ! #variables (desc%nv if 0)
    end subroutine mad_tpsa_poisson

    subroutine mad_tpsa_taylor(tpsa_a,n,coef,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      integer(c_ssz_t), value, intent(in) :: n    ! vector length
      real(c_num_t), intent(in) :: coef(*)        ! vector of taylor coefs
    end subroutine mad_tpsa_taylor

    subroutine mad_tpsa_acc(tpsa_a,v,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! src and dst
      real(c_num_t), value, intent(in) :: v       ! r = r + v*a (r not reset!)
    end subroutine mad_tpsa_acc

    subroutine mad_tpsa_scl(tpsa_a,v,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_num_t), value, intent(in) :: v       ! r = v*a
    end subroutine mad_tpsa_scl

    subroutine mad_tpsa_inv(tpsa_a,v,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_num_t), value, intent(in) :: v       ! r = v/a
    end subroutine mad_tpsa_inv

    subroutine mad_tpsa_invsqrt(tpsa_a,v,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
      real(c_num_t), value, intent(in) :: v       ! r = v/sqrt(a)
    end subroutine mad_tpsa_invsqrt

    subroutine mad_tpsa_sqrt(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sqrt

    subroutine mad_tpsa_exp(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_exp

    subroutine mad_tpsa_log(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_log

    subroutine mad_tpsa_sincos(tpsa_a,tpsa_s,tpsa_c) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_s, tpsa_c        ! dst_sin, dst_cos
    end subroutine mad_tpsa_sincos

    subroutine mad_tpsa_sin(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sin

    subroutine mad_tpsa_cos(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cos

    subroutine mad_tpsa_tan(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_tan

    subroutine mad_tpsa_cot(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cot

    subroutine mad_tpsa_sinc(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinc

    subroutine mad_tpsa_sincosh(tpsa_a,tpsa_s,tpsa_c) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_s, tpsa_c        ! dst_sin, dst_cos
    end subroutine mad_tpsa_sincosh

    subroutine mad_tpsa_sinh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinh

    subroutine mad_tpsa_cosh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_cosh

    subroutine mad_tpsa_tanh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_tanh

    subroutine mad_tpsa_coth(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_coth

    subroutine mad_tpsa_sinhc(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_sinhc

    subroutine mad_tpsa_asin(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_asin

    subroutine mad_tpsa_acos(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acos

    subroutine mad_tpsa_atan(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_atan

    subroutine mad_tpsa_acot(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acot

    subroutine mad_tpsa_asinh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_asinh

    subroutine mad_tpsa_acosh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acosh

    subroutine mad_tpsa_atanh(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_atanh

    subroutine mad_tpsa_acoth(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_acoth

    subroutine mad_tpsa_erf(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_erf

    subroutine mad_tpsa_erfc(tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a    ! src
      type(c_ptr), value :: tpsa_r                ! dst
    end subroutine mad_tpsa_erfc

    ! -- High level functions ---------

    subroutine mad_tpsa_axpb(a,tpsa_x,b,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b    ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x    ! src
      type(c_ptr), value :: tpsa_r                ! dst=a*x+b
    end subroutine mad_tpsa_axpb

    subroutine mad_tpsa_axpbypc(a,tpsa_x,b,tpsa_y,c,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c       ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y  ! src
      type(c_ptr), value :: tpsa_r                      ! dst=a*x+b*y+c
    end subroutine mad_tpsa_axpbypc

    subroutine mad_tpsa_axypb(a,tpsa_x,tpsa_y,b,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b          ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y  ! src
      type(c_ptr), value :: tpsa_r                      ! dst=a*x*y+b
    end subroutine mad_tpsa_axypb

    subroutine mad_tpsa_axypbzpc(a,tpsa_x,tpsa_y,b,tpsa_z,c,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c              ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y, tpsa_z ! src
      type(c_ptr), value :: tpsa_r                             ! dst=a*x*y+b*z+c
    end subroutine mad_tpsa_axypbzpc

    subroutine mad_tpsa_axypbvwpc(a,tpsa_x,tpsa_y,b,tpsa_u,tpsa_v,c,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c            ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y, tpsa_u, tpsa_v ! src
      type(c_ptr), value :: tpsa_r                           ! dst=a*x*y+b*u*v+c
    end subroutine mad_tpsa_axypbvwpc

    subroutine mad_tpsa_ax2pby2pcz2(a,tpsa_x,b,tpsa_y,c,tpsa_z,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c      ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y, tpsa_z ! src
      type(c_ptr), value :: tpsa_r                     ! dst=a*x^2+b*y^2+c*z^2
    end subroutine mad_tpsa_ax2pby2pcz2

    subroutine mad_tpsa_axpsqrtbpcx2(tpsa_x,a,b,c,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c   ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x      ! src
      type(c_ptr), value :: tpsa_r                  ! dst=a*x+sqrt(b+c*x^2)
    end subroutine mad_tpsa_axpsqrtbpcx2

    subroutine mad_tpsa_logaxpsqrtbpcx2(tpsa_x,a,b,c,tpsa_r) bind(C)
      import ; implicit none
      real(c_num_t), value, intent(in) :: a, b, c   ! coefs
      type(c_ptr), value, intent(in) :: tpsa_x      ! src
      type(c_ptr), value :: tpsa_r                  ! dst=log(a*x+sqrt(b+c*x^2))
    end subroutine mad_tpsa_logaxpsqrtbpcx2

    subroutine mad_tpsa_logxdy(tpsa_x,tpsa_y,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_x, tpsa_y ! src
      type(c_ptr), value :: tpsa_r                     ! dst=log(x/y)
    end subroutine mad_tpsa_logxdy

    ! -- Maps-like functions ----------

    subroutine mad_tpsa_minv(n,tpsa_a,tpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n      ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)          ! src
      type(c_ptr) :: tpsa_r(*)                      ! dst
    end subroutine mad_tpsa_minv

    subroutine mad_tpsa_pminv(n,tpsa_a,tpsa_r,select) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n      ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)          ! src
      type(c_ptr) :: tpsa_r(*)                      ! dst
      integer(c_idx_t), intent(in) :: select(*)     ! slots to selected
    end subroutine mad_tpsa_pminv

    subroutine mad_tpsa_compose(na,tpsa_a,nb,tpsa_b,tpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb   ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*), tpsa_b(*) ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst[na]
    end subroutine mad_tpsa_compose

    subroutine mad_tpsa_translate(na,tpsa_a,nb,vb,tpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb   ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      real(c_num_t), intent(in) :: vb(*)              ! src
      type(c_ptr) :: tpsa_r(*)                        ! dst[na]
    end subroutine mad_tpsa_translate

    subroutine mad_tpsa_eval(na,tpsa_a,nb,vb,vr) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb   ! vectors lengths
      type(c_ptr), intent(in) :: tpsa_a(*)            ! src
      real(c_num_t), intent(in) :: vb(*)              ! src
      real(c_num_t) :: vr(*)                          ! dst[nb]
    end subroutine mad_tpsa_eval

    ! -- I/O functions ----------------

    subroutine mad_tpsa_print(tpsa,name_,eps_,nohdr_,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa      ! src
      character(c_char), intent(in) :: name_(*)   ! name (i.e. null terminated string)
      real(c_num_t), value, intent(in) :: eps_    ! display precision, e.g. 1d-12
      integer(c_int), value, intent(in) :: nohdr_ ! discard header if not zero
      type(c_ptr), value :: stream_               ! dst=c_null_ptr => stdio
    end subroutine mad_tpsa_print

    function mad_tpsa_scan(stream_) result (tpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: tpsa                         ! tpsa to read
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end function mad_tpsa_scan

    function mad_tpsa_scan_hdr(kind_,stream_) result (desc) bind(C)
      import ; implicit none
      type(c_ptr) :: desc                         ! descriptor from header
      integer(c_int), optional, intent(out) :: kind_! tpsa kind (0 real, 1 complex)
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end function mad_tpsa_scan_hdr

    subroutine mad_tpsa_scan_coef(tpsa,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value :: tpsa                  ! tpsa to read
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end subroutine mad_tpsa_scan_coef

    subroutine mad_tpsa_debug(tpsa,name_,fnam_,line_,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa     ! src
      character(c_char), intent(in) :: name_(*)  ! name (i.e. null terminated string)
      character(c_char), intent(in) :: fnam_(*)  ! filename (i.e. null terminated string)
      integer(c_int), value, intent(in) :: line_ ! line number or 0
      type(c_ptr), value :: stream_              ! dst=c_null_ptr => stdio
    end subroutine mad_tpsa_debug

  end interface

  ! ----------------------------------------------------------------------------
  ! -- Complex GTPSA -----------------------------------------------------------
  ! ----------------------------------------------------------------------------

  interface

    ! -- Constructors (all constructors return a CTPSA)

    function mad_ctpsa_newd(desc,mo) result(newtpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: newtpsa                    ! new tpsa
      type(c_ptr), value, intent(in) :: desc    ! descriptor
      integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
    end function mad_ctpsa_newd

    function mad_ctpsa_new(ctpsa,mo) result(newtpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: newtpsa                    ! new tpsa
      type(c_ptr), value, intent(in) :: ctpsa   ! (reference) tpsa
      integer(c_ord_t), value, intent(in) :: mo ! if mo > d_mo, mo = d_mo
    end function mad_ctpsa_new

    ! -- Destructor -------------------

    subroutine mad_ctpsa_del(ctpsa) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa   ! tpsa to delete
    end subroutine mad_ctpsa_del

    ! -- Introspection ----------------

    function mad_ctpsa_desc(ctpsa) result(desc) bind(C)
      import ; implicit none
      type(c_ptr) :: desc                       ! tpsa descriptor
      type(c_ptr), value, intent(in) :: ctpsa
    end function mad_ctpsa_desc

    function mad_ctpsa_len(ctpsa) result(len) bind(C)
      import ; implicit none
      integer(c_ssz_t) :: len                   ! #monomials in tpsa
      type(c_ptr), value, intent(in) :: ctpsa
    end function mad_ctpsa_len

    function mad_ctpsa_ord(ctpsa) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                   ! tpsa order
      type(c_ptr), value, intent(in) :: ctpsa
    end function mad_ctpsa_ord

    ! mad_ctpsa_ordv not supported by fortran

    function mad_ctpsa_ordn(n,ctpsa) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                   ! max of all tpsas order
      integer(c_ssz_t), value, intent(in) :: n  ! #ctpsa
      type(c_ptr), intent(in) :: ctpsa(*)
    end function mad_ctpsa_ordn

    function mad_ctpsa_isvalid(ctpsa) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                  ! true or false
      type(c_ptr), value, intent(in) :: ctpsa ! sanity check on TPSA integrity
    end function mad_ctpsa_isvalid

    ! -- Initialization ---------------

    subroutine mad_ctpsa_copy(ctpsa,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_copy

    subroutine mad_ctpsa_getord(ctpsa,ctpsa_r,ord) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
      integer(c_ord_t), value, intent(in) :: ord ! order to retrieve
    end subroutine mad_ctpsa_getord

    subroutine mad_ctpsa_cutord(ctpsa,ctpsa_r,ord) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
      integer(c_int), value, intent(in) :: ord   ! cut order: 0..-ord or ord..mo
    end subroutine mad_ctpsa_cutord

    subroutine mad_ctpsa_convert(ctpsa,ctpsa_r,n,t2r_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa    ! src
      type(c_ptr), value :: ctpsa_r              ! dst
      integer(c_ssz_t), value, intent(in) :: n   ! vector length
      integer(c_idx_t), intent(in) :: t2r_(*)    ! vector of index lookup
    end subroutine mad_ctpsa_convert

    subroutine mad_ctpsa_clear(ctpsa) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
    end subroutine mad_ctpsa_clear

    subroutine mad_ctpsa_setvar(ctpsa,v,iv_,scl_) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      complex(c_cnum_t), value, intent(in) :: v, scl_ ! 0th and 1st order values
      integer(c_idx_t), value, intent(in) :: iv_      ! variable index (1st order)
    end subroutine mad_ctpsa_setvar                   ! equiv. to set0 if iv=0

    ! -- Conversion -------------------

    subroutine mad_ctpsa_real(ctpsa,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa   ! src
      type(c_ptr), value :: tpsa_r              ! dst=real(src)
    end subroutine mad_ctpsa_real

    subroutine mad_ctpsa_imag(ctpsa,tpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa   ! src
      type(c_ptr), value :: tpsa_r              ! dst=imag(src)
    end subroutine mad_ctpsa_imag

    subroutine mad_ctpsa_complex(tpsa_re_,tpsa_im_,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_re_, tpsa_im_ ! src
      type(c_ptr), value :: ctpsa_r                  ! dst=(re or 0)+i*(im or 0)
    end subroutine mad_ctpsa_complex

    ! -- Indexing / monomials (return idx_t = -1 if invalid)

    function mad_ctpsa_mono(ctpsa,n,m_,i) result(ord) bind(C)
      import ; implicit none
      integer(c_ord_t) :: ord                  ! monomial order
      type(c_ptr), value, intent(in) :: ctpsa  !
      integer(c_idx_t), value, intent(in) :: i ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: m_(*)                ! monomial to fill (if provided)
    end function mad_ctpsa_mono

    function mad_ctpsa_idxs(ctpsa,n,s) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index
      type(c_ptr), value, intent(in) :: ctpsa  !
      integer(c_ssz_t), value, intent(in) :: n ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)    ! monomial as string "[0-9]*"
    end function mad_ctpsa_idxs

    function mad_ctpsa_idxm(ctpsa,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index
      type(c_ptr), value, intent(in) :: ctpsa  !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t), intent(in) :: m(*)     ! monomial
    end function mad_ctpsa_idxm

    function mad_ctpsa_idxsm(ctpsa,n,m) result(idx) bind(C)
      import ; implicit none
      integer(c_idx_t) :: idx                  ! monomial index
      type(c_ptr), value, intent(in) :: ctpsa  !
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_int), intent(in) :: m(*)       ! sparse monomial (idx,ord)
    end function mad_ctpsa_idxsm

    function mad_ctpsa_cycle(ctpsa,n,m_,i,v_) result(idx) bind(C)
      import ; implicit none                   ! scan for non-zero coefs starting at i
      integer(c_idx_t) :: idx                  ! next index to start searching or -1
      type(c_ptr), value, intent(in) :: ctpsa  !
      integer(c_idx_t), value, intent(in) :: i ! index to start searching
      integer(c_ssz_t), value, intent(in) :: n ! monomial length
      integer(c_ord_t) :: m_(*)                ! monomial to fill (if provided)
      real(c_cnum_t), intent(out) :: v_        ! coeff to fill (if provided)
    end function mad_ctpsa_cycle

    ! -- Getters ----------------------

    function mad_ctpsa_get0(ctpsa) result(val) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: val                  ! value at order 0 (index 0)
      type(c_ptr), value, intent(in) :: ctpsa
    end function mad_ctpsa_get0

    function mad_ctpsa_geti(ctpsa,i) result(val) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: val                  ! value at index
      type(c_ptr), value, intent(in) :: ctpsa
      integer(c_idx_t), value, intent(in) :: i  ! slot index (must be valid)
    end function mad_ctpsa_geti

    function mad_ctpsa_gets(ctpsa,n,s) result(val) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: val                  ! value at string monomial
      type(c_ptr), value, intent(in) :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n  ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)     ! monomial as string "[0-9]*"
    end function mad_ctpsa_gets

    function mad_ctpsa_getm(ctpsa,n,m) result(val) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: val                  ! value at monomial
      type(c_ptr), value, intent(in) :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n  ! monomial length
      integer(c_ord_t), intent(in) :: m(*)      ! monomial
    end function mad_ctpsa_getm

    function mad_ctpsa_getsm(ctpsa,n,m) result(val) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: val                  ! value at sparse monomial
      type(c_ptr), value, intent(in) :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n  ! monomial length
      integer(c_int), intent(in) :: m(*)        ! sparse monomial (idx,ord)
    end function mad_ctpsa_getsm

    subroutine mad_ctpsa_getv(ctpsa,i,n,v) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa
      integer(c_idx_t), value, intent(in) :: i  ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n  ! vector length
      complex(c_cnum_t) :: v(*)                 ! vector to fill
    end subroutine mad_ctpsa_getv

    ! -- Setters ----------------------

    subroutine mad_ctpsa_set0(ctpsa,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      complex(c_cnum_t), value, intent(in) :: a, b ! ct[0] = a*ct[0]+b
    end subroutine mad_ctpsa_set0

    subroutine mad_ctpsa_seti(ctpsa,i,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      integer(c_idx_t), value, intent(in) :: i      ! slot index (must be valid)
      complex(c_cnum_t), value, intent(in) :: a, b  ! ct[i] = a*ct[i]+b
    end subroutine mad_ctpsa_seti

    subroutine mad_ctpsa_sets(ctpsa,n,s,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n    ! string length or 0 (unknown)
      character(c_char), intent(in) :: s(*)       ! monomial as string "[0-9]*"
      complex(c_cnum_t), value, intent(in) :: a, b ! ct[s] = a*ct[s]+b
    end subroutine mad_ctpsa_sets

    subroutine mad_ctpsa_setm(ctpsa,n,m,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n     ! monomial length
      integer(c_ord_t), intent(in) :: m(*)         ! monomial
      complex(c_cnum_t), value, intent(in) :: a, b ! ct[m] = a*ct[m]+b
    end subroutine mad_ctpsa_setm

    subroutine mad_ctpsa_setsm(ctpsa,n,m,a,b) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      integer(c_ssz_t), value, intent(in) :: n     ! monomial length
      integer(c_int), intent(in) :: m(*)           ! sparse monomial (idx,ord)
      complex(c_cnum_t), value, intent(in) :: a, b ! ct[m] = a*ct[m]+b
    end subroutine mad_ctpsa_setsm

    subroutine mad_ctpsa_setv(ctpsa,i,n,v) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa
      integer(c_idx_t), value, intent(in) :: i      ! slot index (must be valid)
      integer(c_ssz_t), value, intent(in) :: n      ! vector length
      complex(c_cnum_t), intent(in) :: v(*) ! vector to copy
    end subroutine mad_ctpsa_setv

    ! -- Operators --------------------

    function mad_ctpsa_equ(ctpsa_a,ctpsa_b,eps_) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                    ! true or false
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b
      real(c_num_t), value, intent(in) :: eps_  ! tolerance during comparison
    end function mad_ctpsa_equ

    subroutine mad_ctpsa_add(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_add

    subroutine mad_ctpsa_sub(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_sub

    subroutine mad_ctpsa_mul(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_mul

    subroutine mad_ctpsa_div(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_div

    subroutine mad_ctpsa_pow(ctpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                       ! dst
    end subroutine mad_ctpsa_pow

    subroutine mad_ctpsa_powi(ctpsa_a,n,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a   ! src
      integer(c_int), value, intent(in) :: n      ! power (integer)
      type(c_ptr), value :: ctpsa_r               ! dst
    end subroutine mad_ctpsa_powi

    subroutine mad_ctpsa_pown(ctpsa_a,v,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a         ! src
      complex(c_cnum_t), value, intent(in) :: v ! power (real)
      type(c_ptr), value :: ctpsa_r                     ! dst
    end subroutine mad_ctpsa_pown

    ! -- Operators with internal real-to-complex conversion

    function mad_ctpsa_equt(ctpsa_a,tpsa_b,eps) result(ret) bind(C)
      import ; implicit none
      logical(c_bool) :: ret                   ! true or false
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b
      real(c_num_t), value, intent(in) :: eps  ! tolerance during comparison
    end function mad_ctpsa_equt

    subroutine mad_ctpsa_addt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_addt

    subroutine mad_ctpsa_subt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_subt

    subroutine mad_ctpsa_tsub(tpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_tsub

    subroutine mad_ctpsa_mult(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_mult

    subroutine mad_ctpsa_divt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_divt

    subroutine mad_ctpsa_tdiv(tpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_tdiv

    subroutine mad_ctpsa_powt(ctpsa_a,tpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_powt

    subroutine mad_ctpsa_tpow(tpsa_a,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, ctpsa_b  ! lhs, rhs
      type(c_ptr), value :: ctpsa_r                      ! dst
    end subroutine mad_ctpsa_tpow

    subroutine mad_ctpsa_poisst(ctpsa_a,tpsa_b,ctpsa_r,nv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a, tpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int), value, intent(in) :: nv        ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_poisst

    subroutine mad_ctpsa_tpoiss(tpsa_a,ctpsa_b,ctpsa_r,nv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: tpsa_a, ctpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int), value, intent(in) :: nv        ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_tpoiss

    ! -- Functions --------------------

    subroutine mad_ctpsa_abs(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_abs

    subroutine mad_ctpsa_arg(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_arg

    subroutine mad_ctpsa_conj(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a  ! src
      type(c_ptr), value :: ctpsa_r              ! dst
    end subroutine mad_ctpsa_conj

    function mad_ctpsa_nrm1(ctpsa_a,ctpsa_b_) result(nrm1) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: nrm1                   ! sum_i |a[i]-b_[i]|
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b_
    end function mad_ctpsa_nrm1

    function mad_ctpsa_nrm2(ctpsa_a,ctpsa_b_) result(nrm2) bind(C)
      import ; implicit none
      complex(c_cnum_t) :: nrm2                    ! sqrt(sum_i (a[i]-b_[i])^2)
      type(c_ptr), value, intent(in) :: ctpsa_a, ctpsa_b_
    end function mad_ctpsa_nrm2

    subroutine mad_ctpsa_deriv(ctpsa_a,ctpsa_r,iv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      integer(c_int), value, intent(in) :: iv      ! variable index (1st order)
    end subroutine mad_ctpsa_deriv

    subroutine mad_ctpsa_derivm(ctpsa_a,ctpsa_r,n,m) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a   ! src
      type(c_ptr), value :: ctpsa_r               ! dst
      integer(c_ssz_t), value, intent(in) :: n    ! monomial length
      integer(c_ord_t), intent(in) :: m(*)        ! monomial
    end subroutine mad_ctpsa_derivm

    subroutine mad_ctpsa_poisson(ctpsa_a,ctpsa_b,ctpsa_r,nv) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a,ctpsa_b ! src
      type(c_ptr), value :: ctpsa_r                     ! dst
      integer(c_int), value, intent(in) :: nv        ! #variables (desc%nv if 0)
    end subroutine mad_ctpsa_poisson

    subroutine mad_ctpsa_taylor(ctpsa_a,n,coef,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      integer(c_ssz_t), value, intent(in) :: n     ! vector length
      complex(c_cnum_t), intent(in) :: coef(*)     ! vector of taylor coefs
    end subroutine mad_ctpsa_taylor

    subroutine mad_ctpsa_acc(ctpsa_a,v,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! src and dst
      complex(c_cnum_t), value, intent(in) :: v    ! r = r+v*a
    end subroutine mad_ctpsa_acc

    subroutine mad_ctpsa_scl(ctpsa_a,v,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      complex(c_cnum_t), value, intent(in) :: v    ! r = v*a
    end subroutine mad_ctpsa_scl

    subroutine mad_ctpsa_inv(ctpsa_a,v,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      complex(c_cnum_t), value, intent(in) :: v    ! r = v/a
    end subroutine mad_ctpsa_inv

    subroutine mad_ctpsa_invsqrt(ctpsa_a,v,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
      complex(c_cnum_t), value, intent(in) :: v    ! r = v/sqrt(a)
    end subroutine mad_ctpsa_invsqrt

    subroutine mad_ctpsa_sqrt(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sqrt

    subroutine mad_ctpsa_exp(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_exp

    subroutine mad_ctpsa_log(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_log

    subroutine mad_ctpsa_sincos(ctpsa_a,ctpsa_s,ctpsa_c) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a     ! src
      type(c_ptr), value :: ctpsa_s, ctpsa_c        ! dst_sin, dst_cos
    end subroutine mad_ctpsa_sincos

    subroutine mad_ctpsa_sin(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a     ! src
      type(c_ptr), value :: ctpsa_r                 ! dst
    end subroutine mad_ctpsa_sin

    subroutine mad_ctpsa_cos(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cos

    subroutine mad_ctpsa_tan(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_tan

    subroutine mad_ctpsa_cot(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cot

    subroutine mad_ctpsa_sinc(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinc

    subroutine mad_ctpsa_sincosh(ctpsa_a,ctpsa_s,ctpsa_c) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_s, ctpsa_c       ! dst_sin, dst_cos
    end subroutine mad_ctpsa_sincosh

    subroutine mad_ctpsa_sinh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinh

    subroutine mad_ctpsa_cosh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_cosh

    subroutine mad_ctpsa_tanh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_tanh

    subroutine mad_ctpsa_coth(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_coth

    subroutine mad_ctpsa_sinhc(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_sinhc

    subroutine mad_ctpsa_asin(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_asin

    subroutine mad_ctpsa_acos(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acos

    subroutine mad_ctpsa_atan(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_atan

    subroutine mad_ctpsa_acot(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acot

    subroutine mad_ctpsa_asinh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_asinh

    subroutine mad_ctpsa_acosh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acosh

    subroutine mad_ctpsa_atanh(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_atanh

    subroutine mad_ctpsa_acoth(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_acoth

    subroutine mad_ctpsa_erf(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_erf

    subroutine mad_ctpsa_erfc(ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_a    ! src
      type(c_ptr), value :: ctpsa_r                ! dst
    end subroutine mad_ctpsa_erfc

    ! -- High level functions ---------

    subroutine mad_ctpsa_axpb(a,ctpsa_x,b,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x    ! src
      type(c_ptr), value :: ctpsa_r                ! dst=a*x+b
    end subroutine mad_ctpsa_axpb

    subroutine mad_ctpsa_axpbypc(a,ctpsa_x,b,ctpsa_y,c,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c    ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=a*x+b*y+c
    end subroutine mad_ctpsa_axpbypc

    subroutine mad_ctpsa_axypb(a,ctpsa_x,ctpsa_y,b,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b       ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=a*x*y+b
    end subroutine mad_ctpsa_axypb

    subroutine mad_ctpsa_axypbzpc(a,ctpsa_x,ctpsa_y,b,ctpsa_z,c,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c             ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y, ctpsa_z ! src
      type(c_ptr), value :: ctpsa_r                            ! dst=a*x*y+b*z+c
    end subroutine mad_ctpsa_axypbzpc

    subroutine mad_ctpsa_axypbvwpc(a,ctpsa_x,ctpsa_y,b,ctpsa_u,ctpsa_v,c,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c        ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y, ctpsa_u, ctpsa_v ! src
      type(c_ptr), value :: ctpsa_r                          ! dst=a*x*y+b*u*v+c
    end subroutine mad_ctpsa_axypbvwpc

    subroutine mad_ctpsa_ax2pby2pcz2(a,ctpsa_x,b,ctpsa_y,c,ctpsa_z,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c    ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y, ctpsa_z ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=a*x^2+b*y^2+c*z^2
    end subroutine mad_ctpsa_ax2pby2pcz2

    subroutine mad_ctpsa_axpsqrtbpcx2(ctpsa_x,a,b,c,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x       ! src
      type(c_ptr), value :: ctpsa_r                   ! dst=a*x+sqrt(b+c*x^2)
    end subroutine mad_ctpsa_axpsqrtbpcx2

    subroutine mad_ctpsa_logaxpsqrtbpcx2(ctpsa_x,a,b,c,ctpsa_r) bind(C)
      import ; implicit none
      complex(c_cnum_t), value, intent(in) :: a, b, c ! coefs
      type(c_ptr), value, intent(in) :: ctpsa_x       ! src
      type(c_ptr), value :: ctpsa_r                   ! dst=log(a*x+sqrt(b+c*x^2))
    end subroutine mad_ctpsa_logaxpsqrtbpcx2

    subroutine mad_ctpsa_logxdy(ctpsa_x,ctpsa_y,ctpsa_r) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa_x, ctpsa_y ! src
      type(c_ptr), value :: ctpsa_r                      ! dst=log(x/y)
    end subroutine mad_ctpsa_logxdy

    ! -- Maps-like functions ----------

    subroutine mad_ctpsa_minv(n,ctpsa_a,ctpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n        ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)           ! src
      type(c_ptr) :: ctpsa_r(*)                       ! dst
    end subroutine mad_ctpsa_minv

    subroutine mad_ctpsa_pminv(n,ctpsa_a,ctpsa_r,select) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: n        ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)           ! src
      type(c_ptr) :: ctpsa_r(*)                       ! dst
      integer(c_ssz_t), intent(in) :: select(*)       ! slots to selected
    end subroutine mad_ctpsa_pminv

    subroutine mad_ctpsa_compose(na,ctpsa_a,nb,ctpsa_b,ctpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb     ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*), ctpsa_b(*) ! src
      type(c_ptr) :: ctpsa_r(*)                         ! dst[na]
    end subroutine mad_ctpsa_compose

    subroutine mad_ctpsa_translate(na,ctpsa_a,nb,vb,ctpsa_r) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb  ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)          ! src
      complex(c_cnum_t), intent(in) :: vb(*)         ! src
      type(c_ptr) :: ctpsa_r(*)                      ! dst[na]
    end subroutine mad_ctpsa_translate

    subroutine mad_ctpsa_eval(na,ctpsa_a,nb,vb,vr) bind(C)
      import ; implicit none
      integer(c_ssz_t), value, intent(in) :: na, nb  ! vectors lengths
      type(c_ptr), intent(in) :: ctpsa_a(*)          ! src
      complex(c_cnum_t), intent(in) :: vb(*)         ! src
      complex(c_cnum_t) :: vr(*)                     ! dst[nb]
    end subroutine mad_ctpsa_eval

    ! -- I/O functions ----------------

    subroutine mad_ctpsa_print(ctpsa,name_,eps_,nohdr_,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa     ! src
      character(c_char), intent(in) :: name_(*)   ! name (i.e. null terminated string)
      real(c_num_t), value, intent(in) :: eps_    ! display precision, e.g. 1d-12
      integer(c_int), value, intent(in) :: nohdr_ ! discard header if not zero
      type(c_ptr), value :: stream_               ! dst=c_null_ptr => stdio
    end subroutine mad_ctpsa_print

    function mad_ctpsa_scan(stream_) result(ctpsa) bind(C)
      import ; implicit none
      type(c_ptr) :: ctpsa                        ! tpsa to read
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end function mad_ctpsa_scan

    function mad_ctpsa_scan_hdr(kind_,stream_) result(desc) bind(C)
      import ; implicit none
      type(c_ptr) :: desc                         ! descriptor from header
      integer(c_int), optional, intent(out) :: kind_! tpsa kind (0 real, 1 complex)
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end function mad_ctpsa_scan_hdr

    subroutine mad_ctpsa_scan_coef(ctpsa,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value :: ctpsa                 ! tpsa to read
      type(c_ptr), value, intent(in) :: stream_   ! src=c_null_ptr => stdin
    end subroutine mad_ctpsa_scan_coef

    subroutine mad_ctpsa_debug(ctpsa,name_,fnam_,line_,stream_) bind(C)
      import ; implicit none
      type(c_ptr), value, intent(in) :: ctpsa    ! src
      character(c_char), intent(in) :: name_(*)  ! name (i.e. null terminated string)
      character(c_char), intent(in) :: fnam_(*)  ! filename (i.e. null terminated string)
      integer(c_int), value, intent(in) :: line_ ! line number or 0
      type(c_ptr), value :: stream_              ! dst=c_null_ptr => stdio
    end subroutine mad_ctpsa_debug

  end interface

  ! ----------------------------------------------------------------------------
  ! -- C FILEs helpers ---------------------------------------------------------
  ! ----------------------------------------------------------------------------

  interface

    function mad_cio_fopen(path,mode) result(stream) bind(C, name="fopen")
      import ; implicit none
      type(c_ptr) :: stream                             ! see fopen manual
      character(c_char), intent(in) :: path(*), mode(*) ! null terminated strings
    end function mad_cio_fopen

    function mad_cio_fclose(stream) result(ferr) bind(C, name="fclose")
      import ; implicit none
      integer(c_int) :: ferr                            ! file error code or 0
      type(c_ptr), value :: stream                      ! stream to close
    end function mad_cio_fclose

    function mad_cio_fflush(stream) result(ferr) bind(C, name="fflush")
      import ; implicit none
      integer(c_int) :: ferr                            ! file error code or 0
      type(c_ptr), value :: stream                      ! stream to fliush
    end function mad_cio_fflush

    subroutine mad_cio_rewind(stream) bind(C, name="rewind")
      import ; implicit none
      type(c_ptr), value :: stream                      ! stream to rewind
    end subroutine mad_cio_rewind

  end interface

end module GTPSA
