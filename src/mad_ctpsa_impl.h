#ifndef MAD_CTPSA_PRIV_H
#define MAD_CTPSA_PRIV_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Complex Truncated Power Series Algebra module implementation (private)
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
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

#include <tgmath.h>
#include <complex.h>

#include "mad_bit.h"
#include "mad_ctpsa.h"

// --- types ------------------------------------------------------------------o

struct ctpsa { // warning: must be identical to LuaJIT def (see mad_cmad.mad)
  const desc_t *d;
  ord_t lo, hi, mo; // lowest/highest used ord, max ord (allocated)
  bit_t nz;         // warning: must be identical to tpsa up to nz
  cnum_t coef[];
};

// --- macros -----------------------------------------------------------------o

#ifndef MAD_TPSA_NOHELPER

#define T           ctpsa_t
#define NUM         cnum_t
#define FUN(name)   MKNAME(mad_ctpsa_,name)
#define PFX(name)   MKNAME(c,name)
#define VAL(num)    creal(num), cimag(num)
#define FMT         "%+6.4lE%+6.4lEi"
#define SELECT(R,C) C

#define CNUM(a) (* (cnum_t*) & (num_t[2]) { MKNAME(a,_re), MKNAME(a,_im) })

#endif

// --- helpers ----------------------------------------------------------------o

static inline ctpsa_t*
mad_ctpsa_reset0 (ctpsa_t *t)
{
  t->lo = t->mo;
  t->hi = t->nz = 0;
  t->coef[0] = 0;
  return t;
}

static inline ctpsa_t*
mad_ctpsa_copy0 (const ctpsa_t *t, ctpsa_t *r)
{
  if (t != r) {
    r->lo = t->lo;
    r->hi = MIN3(t->hi, r->mo, t->d->to);
    r->nz = mad_bit_hcut(t->nz, r->hi);
    if (r->lo) r->coef[0] = 0;
  }
  return r;
}

static inline ctpsa_t*
mad_ctpsa_update0 (ctpsa_t *t)
{
  idx_t *pi = t->d->ord2idx;
  for (ord_t o = t->lo; o <= t->hi; ++o) {
    if (mad_bit_get(t->nz,o)) {
      idx_t i = pi[o];
      while (i < pi[o+1] && !t->coef[i]) ++i;
      if (i < pi[o+1]) t->nz = mad_bit_clr(t->nz, o);
    }
  }
  return t;
}

static inline ctpsa_t*
mad_ctpsa_gettmp (const ctpsa_t *t, const str_t func)
{
  const D *d = t->d;
  int tid = omp_get_thread_num();
  assert(d->cti[tid] < DESC_MAX_TMP);
  ctpsa_t *tmp = d->ct[ tid*DESC_MAX_TMP + d->cti[tid]++ ];
  TRC_TMPX(printf("GET_TMPX%d[%d]: %p in %s(c)\n",
                  tid, d->cti[tid]-1, (void*)tmp, func));
  tmp->mo = t->mo;
  return tmp; // mad_ctpsa_reset0(tmp);
}

static inline ctpsa_t*
mad_ctpsa_gettmpr (const tpsa_t *t, const str_t func)
{
  return mad_ctpsa_gettmp((const ctpsa_t*)t, func);
}

static inline void
mad_ctpsa_reltmp (ctpsa_t *tmp, const str_t func)
{
  const D *d = tmp->d;
  int tid = omp_get_thread_num();
  TRC_TMPX(printf("REL_TMPX%d[%d]: %p in %s(c)\n",
                  tid, d->cti[tid]-1, (void*)tmp, func));
  assert(d->ct[ tid*DESC_MAX_TMP + d->cti[tid]-1 ] == tmp);
  --d->cti[tid]; //, tmp->mo = d->mo; // ensure stack-like usage of temps
}

// --- end --------------------------------------------------------------------o

#endif // MAD_CTPSA_PRIV_H
