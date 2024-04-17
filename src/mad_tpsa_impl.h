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
  const desc_t *d;      // ptr to tpsa descriptor
  ord_t lo, hi, mo, ao; // lowest/highest used ord, max ord, allocated ord
  int32_t     uid;      // special user field for external use (and padding)
  char  nam[NAMSZ];     // tpsa name (max 15 chars)
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
#define TPSA_SCAN(  ...)   MKNAME(TPSA_SCAN_ ,NARG(__VA_ARGS__))(__VA_ARGS__)
#define TPSA_SCAN_Z(...)   MKNAME(TPSA_SCAN_Z,NARG(__VA_ARGS__))(__VA_ARGS__)
#define TPSA_SCAN_O(...)   MKNAME(TPSA_SCAN_O,NARG(__VA_ARGS__))(__VA_ARGS__)
#define TPSA_SCAN_I(...)   MKNAME(TPSA_SCAN_I,NARG(__VA_ARGS__))(__VA_ARGS__)

#define TPSA_SCAN_1(t)     TPSA_SCAN_3(t,(t)->lo,(t)->hi)
#define TPSA_SCAN_2(t,hi)  TPSA_SCAN_3(t,(t)->lo,     hi)
#define TPSA_SCAN_3(t,lo,hi) \
  const idx_t *o2i = (t)->d->ord2idx; \
  FOR(i,o2i[lo],o2i[(hi)+1]) \

#define TPSA_SCAN_Z1(t)    TPSA_SCAN_Z3(t,(t)->lo,(t)->hi)
#define TPSA_SCAN_Z2(t,hi) TPSA_SCAN_Z3(t,(t)->lo,     hi)
#define TPSA_SCAN_Z3(t,lo,hi) \
  const idx_t *o2i = (t)->d->ord2idx; (void)o2i; \
  for (ord_t o=(lo); o <= (hi); o++) \

#define TPSA_SCAN_I1(t)    TPSA_SCAN_I3(t,(t)->lo,(t)->hi)
#define TPSA_SCAN_I2(t,hi) TPSA_SCAN_I3(t,(t)->lo,     hi)
#define TPSA_SCAN_I3(t,lo,hi) FOR(i,o2i[lo],o2i[(hi)+1])

#define TPSA_SCAN_O1(t)       FOR(i,o2i[o],o2i[o+1])
#define TPSA_SCAN_O2(t,o) \
  const idx_t *o2i = (t)->d->ord2idx; \
  FOR(i,o2i[o],o2i[(o)+1])

// --- functions accessing lo, hi

static inline void // copy TPSA orders but not coefs!
mad_tpsa_copy0 (const tpsa_t *t, tpsa_t *r)
{
  assert(t && r);
  r->lo = t->lo;
  r->hi = MIN(t->hi, r->mo);
  if (r->lo > r->hi) r->lo = 1, r->hi = 0;
}

static inline void // copy TPSA orders but not coefs!
mad_tpsa_copy00 (const tpsa_t *a, const tpsa_t *b, tpsa_t *r)
{
  assert(a && b && r);
  ord_t hi = MAX(a->hi, b->hi);
  r->lo = MIN(a->lo, b->lo);
  r->hi = MIN(hi, r->mo);
  if (r->lo > r->hi) r->lo = 1, r->hi = 0;
}

static inline void // print TPSA header (for debug)
mad_tpsa_print0 (const tpsa_t *a, str_t nam_)
{
  assert(a && a->d);
  printf("'%s' { lo=%d hi=%d mo=%d ao=%d uid=%d did=%d }\n",
         nam_?nam_:"?", a->lo, a->hi, a->mo, a->ao, a->uid, a->d->id);
}

// --- functions accessing lo, hi, coef[0]

static inline tpsa_t* // reset TPSA
mad_tpsa_reset0 (tpsa_t *t)
{
  assert(t);
  t->lo = 1, t->hi = 0, t->coef[0] = 0;
  return t;
}

// --- functions accessing coef[o]

static inline void // clear TPSA order but doesn't adjust lo,hi
mad_tpsa_clear0 (tpsa_t *t, ord_t lo, ord_t hi)
{
  assert(t);
  const idx_t *o2i = t->d->ord2idx;
  FOR(i,o2i[lo],o2i[hi+1]) t->coef[i] = 0;
}

static inline idx_t // return index of first non-zero coef in [lo,hi] or -1
mad_tpsa_nzero0 (const tpsa_t *t, ord_t lo, ord_t hi)
{
  assert(t);
  if (lo > hi) return -1;
  const idx_t *o2i = t->d->ord2idx;
  idx_t i = o2i[lo], ni = o2i[hi+1]-1;
  num_t c = t->coef[ni]; ((tpsa_t*)t)->coef[ni] = 1; // set stopper
  while (!t->coef[i]) ++i;
  ((tpsa_t*)t)->coef[ni] = c;                        // restore value
  return i == ni && !c ? -1 : i;
}

static inline idx_t // return index of first non-zero coef in [lo,hi] or -1
mad_tpsa_nzero0r (const tpsa_t *t, ord_t lo, ord_t hi)
{
  assert(t);
  if (lo > hi) return -1;
  const idx_t *o2i = t->d->ord2idx;
  idx_t ni = o2i[lo], i = o2i[hi+1]-1;
  num_t c = t->coef[ni]; ((tpsa_t*)t)->coef[ni] = 1; // set stopper
  while (!t->coef[i]) --i;
  ((tpsa_t*)t)->coef[ni] = c;                        // restore value
  return i == ni && !c ? -1 : i;
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
