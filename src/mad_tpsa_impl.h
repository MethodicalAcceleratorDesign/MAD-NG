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

struct tpsa_ {  // warning: must be identical to LuaJIT def (see mad_cmad.mad)
  const desc_t *d;  // ptr to tpsa descriptor
  int32_t     uid;  // special user field for external use (and padding)
  ord_t mo, lo, hi; // max ord (allocated), lowest/highest used ord
  bit_t nz;         // zero/non-zero homogeneous polynomials
  char  nam[NAMSZ]; // tpsa name
  num_t coef[]; // warning: must be identical to ctpsa up to coef excluded
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

// loop to scan non-zero homogeneous polynomials (doesn't include coef[0]!)
#define TPSA_SCAN(  ...)      MKNAME(TPSA_SCAN_ ,NARG(__VA_ARGS__))(__VA_ARGS__)
#define TPSA_SCAN_Z(...)      MKNAME(TPSA_SCAN_Z,NARG(__VA_ARGS__))(__VA_ARGS__)
#define TPSA_SCAN_O(...)      MKNAME(TPSA_SCAN_O,NARG(__VA_ARGS__))(__VA_ARGS__)

#define TPSA_SCAN_1(...)      TPSA_SCAN_Z(__VA_ARGS__) TPSA_SCAN_O1()
#define TPSA_SCAN_2(...)      TPSA_SCAN_Z(__VA_ARGS__) TPSA_SCAN_O1()
#define TPSA_SCAN_3(...)      TPSA_SCAN_Z(__VA_ARGS__) TPSA_SCAN_O1()
#define TPSA_SCAN_4(...)      TPSA_SCAN_Z(__VA_ARGS__) TPSA_SCAN_O1()

#define TPSA_SCAN_Z1(t)       TPSA_SCAN_Z4(t,(t)->lo,(t)->hi,(t)->nz)
#define TPSA_SCAN_Z2(t,hi)    TPSA_SCAN_Z4(t,(t)->lo,     hi,(t)->nz)
#define TPSA_SCAN_Z3(t,lo,hi) TPSA_SCAN_Z4(t,     lo,     hi,(t)->nz)
#define TPSA_SCAN_Z4(t,lo,hi,nz) \
  const idx_t *o2i = (t)->d->ord2idx; (void)o2i; \
  const ord_t hi_z = MIN((hi),(t)->d->to); \
  const bit_t nz_z = mad_bit_hcut(nz,hi_z); \
  for (ord_t o = (lo); o <= hi_z; o++) \
    if (mad_bit_tst(nz_z,o))

#define TPSA_SCAN_O1(t) FOR(i,o2i[o],o2i[o+1])
#define TPSA_SCAN_O2(t,o) \
  const idx_t *o2i = (t)->d->ord2idx; \
  FOR(i,o2i[o],o2i[(o)+1])

static inline tpsa_t* // trunc TPSA orders to d->to
mad_tpsa_trunc0 (tpsa_t *t)
{
  assert(t);
  t->hi = MIN(t->hi, t->d->to);
  t->nz = mad_bit_hcut(t->nz, t->hi);
  return t;
}

static inline tpsa_t* // copy TPSA orders but not coefs!
mad_tpsa_copy0 (const tpsa_t *t, tpsa_t *r)
{
  assert(t && r);
  r->lo = MIN(t->lo, r->mo); if (!r->lo) r->lo = 1;
  r->hi = MIN(t->hi, r->mo, r->d->to);
  r->nz = mad_bit_hcut(t->nz, r->hi);
  return r;
}

static inline tpsa_t* // copy TPSA orders but not coefs!
mad_tpsa_copy00 (const tpsa_t *a, const tpsa_t *b, tpsa_t *r)
{
  assert(a && b && r);
  ord_t hi = MAX(a->hi, b->hi);
  r->lo = MIN(a->lo, b->lo, r->mo); if (!r->lo) r->lo = 1;
  r->hi = MIN(hi, r->mo, r->d->to);
  r->nz = mad_bit_add(a->nz, b->nz, r->hi);
  return r;
}

// --- functions accessing coef[0]

static inline log_t // check if TPSA is nul
mad_tpsa_isnul0 (const tpsa_t *t)
{
  assert(t);
  return !(t->nz || t->coef[0]);
}

static inline tpsa_t* // reset TPSA
mad_tpsa_reset0 (tpsa_t *t)
{
  assert(t);
  t->lo = 1, t->hi = 0, t->nz = 0, t->coef[0] = 0;
  return t;
}

static inline tpsa_t* // adjust TPSA orders lo, hi to nz
mad_tpsa_adjust0 (tpsa_t *t)
{
  assert(t);
  if (mad_tpsa_isnul0(t)) return mad_tpsa_reset0(t);
  t->lo = mad_bit_lowest (t->nz);
  t->hi = mad_bit_highest(t->nz);
  return t;
}

// --- functions accessing coef[o]

static inline tpsa_t* // clear TPSA order but doesn't adjust lo,hi
mad_tpsa_clear0 (tpsa_t *t, ord_t o)
{
  assert(t);
  TPSA_SCAN_O(t,o) t->coef[i] = 0;
  return t;
}

static inline tpsa_t* // update TPSA order from coefs but doesn't adjust lo,hi
mad_tpsa_update0 (tpsa_t *t, ord_t o)
{
  assert(t);
  const idx_t *o2i = t->d->ord2idx;
  idx_t i = o2i[o], ni = o2i[o+1]-1;
  num_t c = t->coef[ni]; t->coef[ni] = 1; // set stopper
  while (!t->coef[i]) ++i;
  if (i == ni && !c) t->nz = mad_bit_clr(t->nz,o);
  t->coef[ni] = c; // restore value
  return t;
}

static inline tpsa_t* // round TPSA coefs with magnitude below eps to zero
mad_tpsa_stabilize0 (tpsa_t *t, ord_t o, num_t eps)
{
  assert(t);
  log_t nz = 0;
  TPSA_SCAN_O(t,o)
    if (fabs(t->coef[i]) < eps) t->coef[i] = 0;
    else nz |= !!t->coef[i];
  if (!nz) t->nz = mad_bit_clr(t->nz,o);
  return t;
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
  tmp->mo = t->mo;
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
