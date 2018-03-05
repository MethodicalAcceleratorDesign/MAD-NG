/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map composition module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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
check_same_desc(int sa, const T *ma[sa])
{
  assert(ma);
  for (int i = 1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_compose(int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[])
{
  assert(ma && mb && mc);
  ensure(sa && sb && sc, "invalid map sizes (zero sized map detected)");
  ensure(sa == sc, "incompatibles GTPSA (number of map variables differ)");
  ensure(sb == ma[0]->d->nmv, "incompatibles GTPSA (number of map variables differ)");
  check_same_desc(sa,ma);
  check_same_desc(sb,mb);
  check_same_desc(sc,(const T**)mc);
  ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

#ifdef _OPENMP
#include "mad_tpsa_comp_p.tc"
#endif
#include "mad_tpsa_comp_s.tc"

// --- public -----------------------------------------------------------------o

void
FUN(compose) (int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[])
{
  check_compose(sa, ma, sb, mb, sc, mc);

  #ifdef _OPENMP
  ord_t highest = 0;
  for (int i = 0; i < sa; ++i)
    if (ma[i]->hi > highest) highest = ma[i]->hi;

  if (highest >= 6)
    compose_parallel(sa,ma,mb,mc);
  else
  #endif // _OPENMP

  compose_serial(sa,ma,mb,mc);
}

// --- end --------------------------------------------------------------------o
