/*
 o-----------------------------------------------------------------------------o
 |
 | CTPSA and TPSA conversion module implementation
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

#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#define   MAD_TPSA_NOHELPER
#include "mad_tpsa_impl.h"
#include "mad_ctpsa_impl.h"
#undef    MAD_TPSA_NOHELPER

// --- conversion -------------------------------------------------------------o

void
mad_ctpsa_real (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  const D* d = t->d;
  if (d->to < t->lo) { mad_tpsa_reset0(c); return; }

  c->lo = t->lo;
  c->hi = MIN3(t->hi, c->mo, d->to);
  c->nz = mad_bit_hcut(t->nz,c->hi);
  c->coef[0] = 0;

  idx_t *pi = d->ord2idx;
  for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = creal(t->coef[i]);
}

void
mad_ctpsa_imag (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  const D* d = t->d;
  if (d->to < t->lo) { mad_tpsa_reset0(c); return; }

  c->lo = t->lo;
  c->hi = MIN3(t->hi, c->mo, d->to);
  c->nz = mad_bit_hcut(t->nz,c->hi);
  c->coef[0] = 0;

  idx_t *pi = d->ord2idx;
  for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = cimag(t->coef[i]);
}

void
mad_ctpsa_complex (const tpsa_t *re_, const tpsa_t *im_, ctpsa_t *c)
{
  assert((re_ || im_) && c);
  const tpsa_t *re = re_ ? re_ : im_;
  const tpsa_t *im = im_ ? im_ : re_;
  ensure(re->d == c->d && im->d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t lo = MIN(re->lo, im->lo),
        hi = MAX(re->hi, im->hi);

  const D *d = c->d;
  if (d->to < lo) { mad_ctpsa_reset0(c); return; }

  c->lo = lo;
  c->hi = MIN3(hi, c->mo, d->to);
  c->nz = mad_bit_hcut(re->nz|im->nz,c->hi);
  c->coef[0] = 0;

  idx_t *pi = d->ord2idx;
  switch(!!re_ + 2*!!im_) {
  case 1: // re_ && !im_
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = re->coef[i];
    break;

  case 2: // !re_ && im_
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = im->coef[i]*I;
    break;

  case 3: // re_ && im_
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i) {
      c->coef[i] = 0;
      if (pi[re->lo] <= i && i < pi[re->hi+1]) c->coef[i] += re->coef[i];
      if (pi[im->lo] <= i && i < pi[im->hi+1]) c->coef[i] += im->coef[i]*I;
    }
  }
}

// --- operators with internal conversion -------------------------------------o

log_t mad_ctpsa_equt (const ctpsa_t *a, const tpsa_t *b, num_t tol)
{
  assert(a && b);
  ensure(a->d == b->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  log_t res = mad_ctpsa_equ(a, t, tol);
  REL_TMPC(t);
  return res;
}

void mad_ctpsa_addt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_add(a, t, c);
  REL_TMPC(t);
}

void mad_ctpsa_subt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_sub(a, t, c);
  REL_TMPC(t);
}

void mad_ctpsa_tsub (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(a);
  mad_ctpsa_complex(a, NULL, t);
  mad_ctpsa_sub(t, b, c);
  REL_TMPC(t);
}

void mad_ctpsa_mult (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_mul(a, t, c);
  REL_TMPC(t);
}

void mad_ctpsa_divt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_div(a, t, c);
  REL_TMPC(t);
}

void mad_ctpsa_tdiv (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(a);
  mad_ctpsa_complex(a, NULL, t);
  mad_ctpsa_div(t, b, c);
  REL_TMPC(t);
}

void mad_ctpsa_poisst(const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c, int n)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(b);
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_poisson(a, t, c, n);
  REL_TMPC(t);
}

void mad_ctpsa_tpoiss(const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int n)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPR(a);
  mad_ctpsa_complex(a, NULL, t);
  mad_ctpsa_poisson(t, b, c, n);
  REL_TMPC(t);
}

// --- end --------------------------------------------------------------------o
