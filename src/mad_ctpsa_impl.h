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

struct ctpsa_ { // warning: must be identical to LuaJIT def (see mad_cmad.mad)
  const desc_t *d;  // ptr to ctpsa descriptor
  int32_t     uid;  // special user field for external use (and padding)
  ord_t mo, lo, hi; // max ord (allocated), lowest/highest used ord
  bit_t nz;         // zero/non-zero homogeneous polynomials
  char  nam[NAMSZ]; // tpsa name
  cpx_t coef[]; // warning: must be identical to tpsa up to coef excluded
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

static inline ctpsa_t* // trunc TPSA order to d->to, don't use coef
mad_ctpsa_trunc0 (ctpsa_t *t)
{ return (ctpsa_t*) mad_tpsa_trunc0((tpsa_t*)t); }

static inline ctpsa_t* // adjust TPSA orders lo,hi to nz
mad_ctpsa_adjust0 (ctpsa_t *t)
{ return (ctpsa_t*) mad_tpsa_adjust0((tpsa_t*)t); }

static inline ctpsa_t* // copy TPSA orders, don't use coef.
mad_ctpsa_copy0 (const ctpsa_t *t, ctpsa_t *r)
{ return (ctpsa_t*) mad_tpsa_copy0((const tpsa_t*)t, (tpsa_t*)r); }

static inline ctpsa_t* // copy TPSA orders, don't use coef.
mad_ctpsa_copy00 (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *r)
{ return (ctpsa_t*) mad_tpsa_copy00((const tpsa_t*)a, (const tpsa_t*)b, (tpsa_t*)r); }

// --- functions accessing coef[0]

static inline log_t // check if TPSA is nul
mad_ctpsa_isnul0 (const ctpsa_t *t)
{
  assert(t);
  return !(t->nz || t->coef[0]);
}

static inline ctpsa_t* // reset TPSA
mad_ctpsa_reset0 (ctpsa_t *t)
{
  assert(t);
  t->lo = 1, t->hi = 0, t->nz = 0, t->coef[0] = 0;
  return t;
}

// --- functions accessing coef[o]

static inline ctpsa_t* // clear TPSA order but doesn't adjust lo,hi
mad_ctpsa_clear0 (ctpsa_t *t, ord_t o)
{
  assert(t);
  TPSA_SCAN_O(t,o) t->coef[i] = 0;
  return t;
}

static inline ctpsa_t* // update TPSA order from coefs but doesn't adjust lo,hi
mad_ctpsa_update0 (ctpsa_t *t, ord_t o)
{
  assert(t);
  const idx_t *o2i = t->d->ord2idx;
  idx_t i = o2i[o], ni = o2i[o+1]-1;
  cpx_t c = t->coef[ni]; t->coef[ni] = 1; // set stopper
  while (!t->coef[i]) ++i;
  if (i == ni && !c) t->nz = mad_bit_clr(t->nz,o);
  t->coef[ni] = c; // restore value
  return t;
}

static inline ctpsa_t* // round TPSA coefs with magnitude below eps to zero
mad_ctpsa_stabilize0 (ctpsa_t *t, ord_t o, num_t eps)
{
  assert(t);
  log_t nz = 0;
  TPSA_SCAN_O(t,o)
    if (fabs(creal(t->coef[i])) < eps &&
        fabs(cimag(t->coef[i])) < eps) t->coef[i] = 0;
    else nz |= !!t->coef[i];
  if (!nz) t->nz = mad_bit_clr(t->nz,o);
  return t;
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
