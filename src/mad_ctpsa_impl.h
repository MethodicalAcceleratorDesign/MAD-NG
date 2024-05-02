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

#include <assert.h>
#include <tgmath.h>
#include <complex.h>

#include "mad_bit.h"
#include "mad_ctpsa.h"
#include "mad_tpsa_impl.h"

// --- types ------------------------------------------------------------------o

struct ctpsa_ { // warning: must be identical to LuaJIT def (see mad_gtpsa.mad)
  const desc_t *d;        // ptr to ctpsa descriptor
  ord_t   lo, hi, mo, ao; // lowest/highest used ord, max ord, allocated ord
  int32_t uid;            // special user field for external use (and padding)
  char    nam[NAMSZ];     // tpsa name (max 15 chars)
  cpx_t   coef[]; // warning: must be identical to tpsa up to coef excluded
};

// --- macros -----------------------------------------------------------------o

#ifdef MAD_CTPSA_IMPL

#define T                ctpsa_t
#define NUM              cpx_t
#define FUN(name)        MKNAME(mad_ctpsa_,name)
#define PFX(name)        MKNAME(c,name)
#define VAL(num)         creal(num), cimag(num)
#define VALEPS(num,eps) (fabs(creal(num))<(eps) ? 0 : creal(num)), \
                        (fabs(cimag(num))<(eps) ? 0 : cimag(num))
#define FMT              "%+6.4lE%+6.4lEi"
#define SELECT(R,C)      C

#define CPX(a) (* (cpx_t*) & (num_t[2]) { MKNAME(a,_re), MKNAME(a,_im) })

#endif

// --- helpers ----------------------------------------------------------------o

// --- functions accessing lo, hi

static inline void // copy TPSA orders, don't use coef.
mad_ctpsa_copy0 (const ctpsa_t *t, ctpsa_t *r)
{ mad_tpsa_copy0((const tpsa_t*)t, (tpsa_t*)r); }

static inline void // copy TPSA orders, don't use coef.
mad_ctpsa_copy00 (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *r)
{ mad_tpsa_copy00((const tpsa_t*)a, (const tpsa_t*)b, (tpsa_t*)r); }

static inline void // print TPSA header (for debug)
mad_ctpsa_print0 (const ctpsa_t *t, str_t nam_)
{ mad_tpsa_print0((const tpsa_t*)t, nam_); }

// --- functions accessing lo, hi, coef[0]

static inline ctpsa_t* // reset TPSA
mad_ctpsa_reset0 (ctpsa_t *t)
{
  assert(t);
  t->lo = 1, t->hi = 0, t->coef[0] = 0;
  return t;
}

// --- functions accessing coef[o]

static inline void // clear TPSA order but doesn't adjust lo,hi
mad_ctpsa_clear0 (ctpsa_t *t, ord_t lo, ord_t hi)
{
  assert(t);
  const idx_t *o2i = t->d->ord2idx;
  FOR(i,o2i[lo],o2i[hi+1]) t->coef[i] = 0;
}

static inline idx_t // return index of first non-zero coef in [lo,hi] or -1
mad_ctpsa_nzero0 (const ctpsa_t *t, ord_t lo, ord_t hi, log_t upt)
{
  assert(t);
  if (lo > hi) return -1;
  const idx_t *o2i = t->d->ord2idx;
  idx_t i = o2i[lo], ni = o2i[hi+1]-1;
  cpx_t c = t->coef[ni]; ((ctpsa_t*)t)->coef[ni] = 1; // set stopper
  while (!t->coef[i]) ++i;
  ((ctpsa_t*)t)->coef[ni] = c;                        // restore value
  if (i != ni || c) {
    if (upt) ((ctpsa_t*)t)->lo = t->d->ords[i];
    return i;
  } else {
    if (upt) ((ctpsa_t*)t)->lo = 1, ((ctpsa_t*)t)->hi = 0;
    return -1;
  }
}

static inline idx_t // return index of first non-zero coef in [lo,hi] or -1
mad_ctpsa_nzero0r (const ctpsa_t *t, ord_t lo, ord_t hi, log_t upt)
{
  assert(t);
  const idx_t *o2i = t->d->ord2idx;
  for (ord_t o = hi; o >= lo; o--) {
    idx_t i = o2i[o], ni = o2i[o+1]-1;
    cpx_t c = t->coef[ni]; ((ctpsa_t*)t)->coef[ni] = 1; // set stopper
    while (!t->coef[i]) ++i;
    ((ctpsa_t*)t)->coef[ni] = c;                        // restore value
    if (i != ni || c) {
      if (upt) ((ctpsa_t*)t)->hi = o;
      return i;
    }
  }
  return -1;
}

// --- temporaries ------------------------------------------------------------o

#if DESC_USE_TMP

static inline ctpsa_t*
mad_ctpsa_gettmp (const ctpsa_t *t, const str_t func)
{
  assert(t);
  const desc_t *d = t->d;
  int tid = omp_get_thread_num();
  assert(d->cti[tid] < DESC_MAX_TMP);
  ctpsa_t *tmp = d->ct[ tid*DESC_MAX_TMP + d->cti[tid]++ ];
  TRC_TMPX(printf("GET_TMPX%d[%d]: %p in %s(c)\n",
                  tid, d->cti[tid]-1, (void*)tmp, func));
  tmp->mo = t->mo;
  return mad_ctpsa_reset0(tmp);
}

static inline void
mad_ctpsa_reltmp (ctpsa_t *tmp, const str_t func)
{
  assert(tmp);
  const desc_t *d = tmp->d;
  int tid = omp_get_thread_num();
  TRC_TMPX(printf("REL_TMPX%d[%d]: %p in %s(c)\n",
                  tid, d->cti[tid]-1, (void*)tmp, func));
  assert(d->ct[ tid*DESC_MAX_TMP + d->cti[tid]-1 ] == tmp);
  --d->cti[tid]; //, tmp->mo = d->mo; // ensure stack-like usage of temps
}

static inline ctpsa_t*
mad_ctpsa_gettmpt (const tpsa_t *t, const str_t func)
{
  return mad_ctpsa_gettmp((const ctpsa_t*)t, func);
}

#endif // DESC_USE_TMP

// --- end --------------------------------------------------------------------o

#endif // MAD_CTPSA_PRIV_H
