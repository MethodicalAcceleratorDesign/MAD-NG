/*
 o----------------------------------------------------------------------------o
 |
 | TPSA map composition module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

#include <math.h>
#include <string.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_tpsa.h"

#include "mad_tpsa_impl.h"
#include "mad_desc_impl.h"

#define T struct tpsa
#define D struct tpsa_desc

#undef  ensure
#define ensure(test) mad_ensure(test, MKSTR(test))

// --- LOCAL FUNCTIONS --------------------------------------------------------

static inline void
check_same_desc(int sa, const T *ma[sa])
{
  assert(ma);
  for (int i = 1; i < sa; ++i)
    ensure(ma[i]->desc == ma[i-1]->desc);
}

static inline void
check_compose(int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[])
{
  assert(ma && mb && mc);
  ensure(sa && sb && sc);
  ensure(sa == sc);
  ensure(sb == ma[0]->desc->nmv);
  check_same_desc(sa,ma);
  check_same_desc(sb,mb);
  check_same_desc(sc,(const T**)mc);
  ensure(ma[0]->desc == mb[0]->desc);
  ensure(ma[0]->desc == mc[0]->desc);
}

#ifdef _OPENMP
#include "mad_tpsa_comp_p.tc"
#endif
#include "mad_tpsa_comp_s.tc"

// --- PUBLIC FUNCTIONS -------------------------------------------------------

void
mad_tpsa_compose(int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[])
{
#ifdef TRACE
  printf("tpsa_compose\n");
#endif
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
