#ifndef MAD_TPSA_PRIV_H
#define MAD_TPSA_PRIV_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation (private)
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
#include <assert.h>

#include "mad_bit.h"
#include "mad_tpsa.h"
#include "mad_desc_impl.h"

// --- types ------------------------------------------------------------------o

struct tpsa_ {  // warning: must be identical to LuaJIT def (see mad_gtpsa.mad)
  const desc_t *d;     // ptr to tpsa descriptor
  ord_t   mo, ao;      // max ord, allocated ord
  ssz_t   nc, mc;      // #val, max #val = index of mo+1
  int32_t uid;         // special "user" field for external use (padding)
  char    nam[NAMSZ];  // tpsa name (hidden)
  idx_t  *idx;         // = (idx_t*)(val+mc);
  num_t   val[];       // R:nn+(nn+1)/2, C:nn+(nn+2)/4, no = ords[nn]
};

// --- macros -----------------------------------------------------------------o

#ifndef MAD_CTPSA_IMPL

#define T                tpsa_t
#define NUM              num_t
#define FUN(name)        MKNAME(mad_tpsa_,name)
#define PFX(name)        name
#define VAL(num)         num
#define VALEPS(num,eps) (fabs(num)<(eps) ? 0 : (num))
#define FMT              "%+6.4lE"
#define SELECT(R,C)      R

#endif

// --- helpers ----------------------------------------------------------------o

static inline tpsa_t* // reset TPSA
mad_tpsa_reset0 (tpsa_t *t)
{
  assert(t);
  t->nc = 1, t->val[0] = 0;
  return t;
}

static inline void // print TPSA header (for debug)
mad_tpsa_print0 (const tpsa_t *t, str_t nam_)
{
  assert(t && t->d);
  printf("'%s' { nc=%d mc=%d mo=%d ao=%d uid=%d did=%d }\n",
         nam_?nam_:"?", t->nc, t->mc, t->mo, t->ao, t->uid, t->d->id);
}

// --- temporaries ------------------------------------------------------------o

#if DESC_USE_TMP

static inline tpsa_t*
mad_tpsa_gettmp (const tpsa_t *t, const str_t func)
{
  assert(t);
  const desc_t *d = t->d;
  int tid = omp_get_thread_num();
  assert(d->ti[tid] < DESC_MAX_TMP);
  tpsa_t *tmp = d->t[ tid*DESC_MAX_TMP + d->ti[tid]++ ];
  TRC_TMPX(printf("GET_TMPX%d[%d]: %p in %s()\n",
                  tid, d->ti[tid]-1, (void*)tmp, func));
  tmp->mo = t->mo, tmp->mc = t->mc;
  return mad_tpsa_reset0(tmp);
}

static inline void
mad_tpsa_reltmp (tpsa_t *tmp, const str_t func)
{
  assert(tmp);
  const desc_t *d = tmp->d;
  int tid = omp_get_thread_num();
  TRC_TMPX(printf("REL_TMPX%d[%d]: %p in %s()\n",
                  tid, d->ti[tid]-1, (void*)tmp, func));
  assert(d->t[ tid*DESC_MAX_TMP + d->ti[tid]-1 ] == tmp);
  --d->ti[tid]; // ensure stack-like usage of temps
}

static inline tpsa_t*
mad_tpsa_gettmpt (const ctpsa_t *t, const str_t func)
{
  return mad_tpsa_gettmp((const tpsa_t*)t, func);
}

#endif // DESC_USE_TMP

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_PRIV_H
