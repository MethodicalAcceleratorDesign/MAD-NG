/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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

  desc_t *d = t->d;
  if (d->trunc < t->lo) { mad_tpsa_clear(c); return; }

  c->lo = t->lo;
  c->hi = MIN3(t->hi, c->mo, d->trunc);
  c->nz = mad_bit_trunc(t->nz, c->hi);

  idx_t *pi = d->ord2idx;

  for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = creal(t->coef[i]);
}

void
mad_ctpsa_imag (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  desc_t *d = t->d;
  if (d->trunc < t->lo) { mad_tpsa_clear(c); return; }

  c->lo = t->lo;
  c->hi = MIN3(t->hi, c->mo, d->trunc);
  c->nz = mad_bit_trunc(t->nz, c->hi);

  idx_t *pi = d->ord2idx;

  for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
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

  D *d = c->d;
  if (d->trunc < lo) { mad_ctpsa_clear(c); return; }

  c->lo  = lo;
  c->hi  = MIN3(hi, c->mo, d->trunc);
  c->nz  = mad_bit_trunc(mad_bit_add(re->nz, im->nz), c->hi);

  idx_t *pi = d->ord2idx;

  switch(!!re_ + 2*!!im_) {
  case 1: // re_ && !im_
    for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = re->coef[i];
    break;

  case 2: // !re_ && im_
    for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = im->coef[i]*I;
    break;

  case 3: // re_ && im_
    for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i) {
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
  ctpsa_t *t = a->d->ct[0];
  mad_ctpsa_complex(b, NULL, t);
  return mad_ctpsa_equ(a, t, tol);
}

void mad_ctpsa_addt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[0];
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_add(a, t, c);
}

void mad_ctpsa_subt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[0];
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_sub(a, t, c);
}

void mad_ctpsa_tsub (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[0];
  mad_ctpsa_complex(a, NULL, t);
  mad_ctpsa_sub(t, b, c);
}

void mad_ctpsa_mult (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[1];      // mul.tmp uses ct[0] if a == c
  mad_ctpsa_complex(b, NULL, t);
  mad_ctpsa_mul(a, t, c);
}

void mad_ctpsa_divt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[0];      // div.tmp uses ct[4]
  mad_ctpsa_complex(b, NULL, t); // div.inv uses fix pts ct[1-3]
  mad_ctpsa_div(a, t, c);        // div.mul uses ct[0] if a == c, but not t (!)
}

void mad_ctpsa_tdiv (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[0];      // see divt
  mad_ctpsa_complex(a, NULL, t);
  mad_ctpsa_div(t, b, c);
}

void mad_ctpsa_poisst(const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c, int n)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[1];
  mad_ctpsa_complex(b, NULL, t); // poisson.mul uses ct[0]
  mad_ctpsa_poisson(a, t, c, n);
}

void mad_ctpsa_tpoiss(const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int n)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = c->d->ct[1];
  mad_ctpsa_complex(a, NULL, t); // poisson.mul uses ct[0]
  mad_ctpsa_poisson(t, b, c, n);
}

// --- end --------------------------------------------------------------------o
