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
  const desc_t *d;     // ptr to tpsa descriptor
  ord_t   mo, ao;      // max ord, allocated ord
  ssz_t   nc, mc;      // #val, max #val = index of mo+1
  int32_t uid;         // special "user" field for external use (padding)
  char    nam[NAMSZ];  // tpsa name (hidden)
  idx_t  *idx;         // = (idx_t*)(val+nn);
  cpx_t   val[];       // R:nn+(nn+1)/2, C:nn+(nn+2)/4, no = ords[nn]
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

static inline ctpsa_t* // reset TPSA
mad_ctpsa_reset0 (ctpsa_t *t)
{
  assert(t);
  t->nc = 1, t->val[0] = 0;
  return t;
}

static inline void // print TPSA header (for debug)
mad_ctpsa_print0 (const ctpsa_t *t, str_t nam_)
{
  assert(t && t->d);
  printf("'%s' { nc=%d mc=%d mo=%d ao=%d uid=%d did=%d }\n",
         nam_?nam_:"?", t->nc, t->mc, t->mo, t->ao, t->uid, t->d->id);
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
  tmp->mo = t->mo, tmp->mc = t->mc;
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
