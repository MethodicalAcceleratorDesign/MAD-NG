/*
 o-----------------------------------------------------------------------------o
 |
 | CTPSA map normal form module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
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
#include "mad_ctpsa_impl.h"

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc(ssz_t sa, const T *ma[sa])
{
  assert(ma);
  for (idx_t i = 1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_normal(ssz_t sa, const T *ma[sa], ssz_t sb, T *mb[sb], ssz_t sc, T *mc[sc])
{
  assert(ma && mb && mc);
  ensure(sa>0 && sb>0 && sc>0, "invalid map sizes (zero or negative sized map)");
  ensure(sa == sc, "incompatibles GTPSA (number of map variables differ)");
  ensure(sb == ma[0]->d->nmv, "incompatibles GTPSA (number of map variables differ)");
  check_same_desc(sa,ma);
  check_same_desc(sb,(const T**)mb);
  check_same_desc(sc,(const T**)mc);
  ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

/* TODO
#ifdef _OPENMP
#include "mad_tpsa_nform_p.tc"
#endif
#include "mad_tpsa_nform_s.tc"
*/

// --- public -----------------------------------------------------------------o

// map normal form
// ma = mc mb mc^-1,
// x  =  a  r  a^-1, where r is the normalized map (containing tunes and dampings)

void
FUN(normal) (ssz_t sa, const T *ma[sa], ssz_t sb, T *mb[sb], ssz_t sc, T *mc[sc])
{
  ensure(0, "NYI: high order normal form");
  check_normal(sa, ma, sb, mb, sc, mc);

/* TODO
  #ifdef _OPENMP
  ord_t highest = 0;
  for (idx_t i = 0; i < sa; ++i)
    if (ma[i]->hi > highest) highest = ma[i]->hi;

  if (highest >= 6)
    normal_parallel(sa,ma,mb,mc);
  else
  #endif // _OPENMP

  normal_serial(sa,ma,mb,mc);
*/
}

void
FUN(tnormal)(ssz_t sa, const tpsa_t *ma[sa], ssz_t sb, T *mb[sb], ssz_t sc, T *mc[sc])
{
  ensure(0, "NYI: high order normal form");
  check_normal(sa, (const T**)ma, sb, mb, sc, mc);

  const D *d = mb[0]->d;
  const T *cma[sa];
  T *a;

#undef complex // undefined macro complex from <complex.h>

  for (idx_t i = 0; i < sa; i++) {
    a = FUN(newd)(d, mad_tpsa_same);
    FUN(complex)(ma[i], 0, a);
    cma[i] = a;
  }

  FUN(normal)(sa, cma, sb, mb, sc, mc);

  for (idx_t i = 0; i < sa; i++)
    FUN(del)(cma[i]);
}

// --- end --------------------------------------------------------------------o
