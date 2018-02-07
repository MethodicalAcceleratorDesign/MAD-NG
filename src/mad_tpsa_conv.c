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

void
mad_ctpsa_real (const ctpsa_t *t, tpsa_t *dst)
{
  assert(t && dst);
  ensure(t->d == dst->d, "incompatibles GTPSA (descriptors differ)");

  desc_t *d = t->d;
  if (d->trunc < t->lo) { mad_tpsa_clear(dst); return; }

  dst->lo = t->lo;
  dst->hi = MIN3(t->hi, dst->mo, d->trunc);
  dst->nz = mad_bit_trunc(t->nz, dst->hi);

  for (int i = d->ord2idx[dst->lo]; i < d->ord2idx[dst->hi+1]; ++i)
    dst->coef[i] = creal(t->coef[i]);
}

void
mad_ctpsa_imag (const ctpsa_t *t, tpsa_t *dst)
{
  assert(t && dst);
  ensure(t->d == dst->d, "incompatibles GTPSA (descriptors differ)");

  desc_t *d = t->d;
  if (d->trunc < t->lo) { mad_tpsa_clear(dst); return; }

  dst->lo = t->lo;
  dst->hi = MIN3(t->hi, dst->mo, d->trunc);
  dst->nz = mad_bit_trunc(t->nz, dst->hi);

  for (int i = d->ord2idx[dst->lo]; i < d->ord2idx[dst->hi+1]; ++i)
    dst->coef[i] = cimag(t->coef[i]);
}

void
mad_tpsa_complex (const tpsa_t *re_, const tpsa_t *im_, ctpsa_t *dst)
{
  assert((re_ || im_) && dst);
  const tpsa_t *re = re_ ? re_ : im_;
  const tpsa_t *im = im_ ? im_ : re_;
  ensure(re->d == dst->d && im->d == dst->d, "incompatibles GTPSA (descriptors differ)");

  D *d = dst->d;
  if (d->trunc < MIN(re->lo, im->lo)) { mad_ctpsa_clear(dst); return; }

  ord_t hi = MAX(re->hi, im->hi);
  dst->lo  = MIN(re->lo, im->lo);
  dst->hi  = MIN3(hi, dst->mo, d->trunc);
  dst->nz  = mad_bit_trunc(mad_bit_add(re->nz, im->nz), dst->hi);

  for (int i = d->ord2idx[dst->lo]; i < d->ord2idx[dst->hi+1]; ++i) {
    dst->coef[i] = 0;
    if (re_ && d->ord2idx[re->lo] <= i && i < d->ord2idx[re->hi+1])
      dst->coef[i] += re->coef[i];
    if (im_ && d->ord2idx[im->lo] <= i && i < d->ord2idx[im->hi+1])
      dst->coef[i] += im->coef[i]*I;
  }
}

// --- end --------------------------------------------------------------------o
