/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map composition module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o
*/

#include <math.h>
#include <string.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc (ssz_t sa, const T *ma[sa])
{
  assert(ma);
  for (idx_t i = 1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_compose (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  assert(ma && mb && mc);
  ensure(sa>0 && sb>0, "invalid map sizes (zero or negative sizes)");
  ensure(sb == ma[0]->d->nmv, "incompatibles GTPSA (number of map variables differ)");
  check_same_desc(sa,ma);
  check_same_desc(sb,mb);
  check_same_desc(sa,(const T**)mc);
  ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
print_damap (ssz_t sa, const T *ma[sa])
{
  char nam[3] = "#1";
  for (ssz_t i=0; i < sa; i++, nam[1]++) FUN(print)(ma[i], nam, 0, 0, stdout);
}

#ifdef _OPENMP
#include "mad_tpsa_comp_p.tc"
#endif
#include "mad_tpsa_comp_s.tc"

// --- public -----------------------------------------------------------------o

void
FUN(compose) (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  DBGFUN(->);
  check_compose(sa, ma, sb, mb, mc);

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  for (idx_t ic = 0; ic < sa; ++ic) {
    DBGTPSA(ma[ic]); DBGTPSA(mc[ic]);
    mc_[ic] = FUN(new)(mc[ic], mad_tpsa_same);
  }

  for (idx_t ib = 0; ib < sb; ++ib) {
    DBGTPSA(mb[ib]);
  }

  #ifdef _OPENMP
  ord_t highest = 0;
  for (idx_t i = 0; i < sa; ++i)
    if (ma[i]->hi > highest) highest = ma[i]->hi;

  if (highest >= 6)
    compose_parallel(sa,ma,mb,mc_);
  else
  #endif // _OPENMP
    compose_serial  (sa,ma,mb,mc_);

  // copy back
  for (idx_t ic = 0; ic < sa; ++ic) {
    FUN(copy)(mc_[ic], mc[ic]);
    FUN(del )(mc_[ic]);
    DBGTPSA(mc[ic]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void
FUN(translate) (ssz_t sa, const T *ma[sa], ssz_t sb, const NUM tb[sb], T *mc[sa])
{
  assert(ma && tb); DBGFUN(->);
  ensure(sb>0, "invalid vector sizes (zero or negative sizes)");

  // transform translation vector into damap of order 1
  mad_alloc_tmp(const T*, mb, sb);
  for (idx_t ib = 0; ib < sb; ++ib) {
    T *t = FUN(newd)(ma[0]->d, 1);
    FUN(setvar)(t, tb[ib], ib+1, 0);
    mb[ib] = t;
  }

  FUN(compose)(sa, ma, sb, mb, mc);

  // cleanup
  for (idx_t ib = 0; ib < sb; ++ib)
    FUN(del)(mb[ib]);
  mad_free_tmp(mb);
  DBGFUN(<-);
}

void
FUN(eval) (ssz_t sa, const T *ma[sa], ssz_t sb, const NUM tb[sb], NUM tc[sb])
{
  assert(ma && tb && tc); DBGFUN(->);
  ensure(sa>0 && sb>0, "invalid map/vector sizes (zero or negative sizes)");
  ensure(sb == ma[0]->d->nmv, "incompatibles GTPSA (number of map variables differ)");

  // transform vectors into damap of order 0
  mad_alloc_tmp(const T*, mb, sb);
  mad_alloc_tmp(      T*, mc, sb);
  for (idx_t ib = 0; ib < sb; ++ib) {
    T *t1 = FUN(newd)(ma[0]->d, 0);
    T *t2 = FUN(newd)(ma[0]->d, 0);
    FUN(setvar)(t1, tb[ib], 0, 0);
    FUN(setvar)(t2, tc[ib], 0, 0);
    mb[ib] = t1; mc[ib] = t2;
  }

  compose_serial(sa,ma,mb,mc);

  // cleanup, save result
  for (idx_t ib = 0; ib < sb; ++ib) {
    tc[ib] = mc[ib]->coef[0];
    FUN(del)(mc[ib]);
    FUN(del)(mb[ib]);
  }
  mad_free_tmp(mc);
  mad_free_tmp(mb);
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
