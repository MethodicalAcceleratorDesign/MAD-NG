/*
 o-----------------------------------------------------------------------------o
 |
 | CTPSA and TPSA conversion module implementation
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

#include "mad_mem.h"
#include "mad_desc_impl.h"

#include "mad_tpsa_impl.h"
#define   MAD_TPSA_NOHELPER
#include "mad_ctpsa_impl.h"
#undef    MAD_TPSA_NOHELPER

// --- conversion -------------------------------------------------------------o

void
mad_ctpsa_real (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c); DBGFUN(->);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  const D* d = t->d;

  c->hi = MIN3(t->hi, c->mo, d->to);
  c->nz = mad_bit_hcut(t->nz,c->hi);
  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  c->lo = t->lo;
  c->coef[0] = 0;

  const idx_t *o2i = d->ord2idx;
  for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
    c->coef[i] = creal(t->coef[i]);

  if (TPSA_STRICT_NZ) FUN(update0)(c, c->lo, c->hi);

  DBGFUN(<-);
}

void
mad_ctpsa_imag (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c); DBGFUN(->);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  const D* d = t->d;

  c->hi = MIN3(t->hi, c->mo, d->to);
  c->nz = mad_bit_hcut(t->nz,c->hi);
  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  c->lo = t->lo;
  c->coef[0] = 0;

  const idx_t *o2i = d->ord2idx;
  for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
    c->coef[i] = cimag(t->coef[i]);

  if (TPSA_STRICT_NZ) FUN(update0)(c, c->lo, c->hi);

  DBGFUN(<-);
}

void // special unique case, should be in mad_tpsa_conv.c to use FUN
mad_ctpsa_cplx (const tpsa_t *re_, const tpsa_t *im_, ctpsa_t *c)
{
  assert((re_ || im_) && c); DBGFUN(->);
  const tpsa_t *re = re_ ? re_ : im_;
  const tpsa_t *im = im_ ? im_ : re_;
  ensure(re->d == c->d && im->d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t lo = MIN(re->lo, im->lo),
        hi = MAX(re->hi, im->hi);

  const D *d = c->d;

  c->hi = MIN3(hi, c->mo, d->to);
  c->nz = mad_bit_hcut(re->nz|im->nz,c->hi);
  if (!c->nz) { mad_ctpsa_reset0(c); DBGFUN(<-); return; }

  c->lo = lo;
  c->coef[0] = 0;

  const idx_t *o2i = d->ord2idx;
  switch(!!re_ + 2*!!im_) {
  case 1: // re_ && !im_
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
      c->coef[i] = re_->coef[i];
    break;

  case 2: // !re_ && im_
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
      c->coef[i] = im_->coef[i]*I;
    break;

  case 3: // re_ && im_
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i) {
      c->coef[i] = 0;
      if (o2i[re_->lo] <= i && i < o2i[re_->hi+1]) c->coef[i] += re_->coef[i];
      if (o2i[im_->lo] <= i && i < o2i[im_->hi+1]) c->coef[i] += im_->coef[i]*I;
    }
  }
  DBGFUN(<-);
}

// --- operators with internal conversion -------------------------------------o

void mad_ctpsa_cabs (const ctpsa_t *a, tpsa_t *c)
{
  assert(a && c); DBGFUN(->);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  tpsa_t *re = GET_TMPR(a), *im = GET_TMPR(a);
  mad_ctpsa_real(a, re);
  mad_ctpsa_imag(a, im);
  mad_tpsa_hypot(re, im, c);
  REL_TMPR(re); REL_TMPR(im);
  DBGFUN(<-);
}

void mad_ctpsa_carg (const ctpsa_t *a, tpsa_t *c)
{
  assert(a && c); DBGFUN(->);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  tpsa_t *re = GET_TMPR(a), *im = GET_TMPR(a);
  mad_ctpsa_real(a, re);
  mad_ctpsa_imag(a, im);
  mad_tpsa_atan2(im, re, c);
  REL_TMPR(re); REL_TMPR(im);
  DBGFUN(<-);
}

void mad_ctpsa_rect (const ctpsa_t *a, ctpsa_t *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  tpsa_t *re = GET_TMPR(c), *im = GET_TMPR(c), *t = GET_TMPR(c);

  mad_ctpsa_real (a , re); // rho
  mad_ctpsa_imag (a , im); // arg
  mad_tpsa_sincos(im, im, t );
  mad_tpsa_mul   (re, t , re); // re = rho * cos(arg)
  mad_tpsa_mul   (re, im, im); // im = rho * sin(arg)
  mad_ctpsa_cplx (re, im, c );

  REL_TMPR(re); REL_TMPR(im); REL_TMPR(t);

  DBGTPSA(c); DBGFUN(<-);
}

void mad_ctpsa_polar (const ctpsa_t *a, ctpsa_t *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  tpsa_t *re = GET_TMPR(c), *im = GET_TMPR(c), *t = GET_TMPR(c);

  mad_ctpsa_real(a , re);
  mad_ctpsa_imag(a , im);
  mad_tpsa_hypot(im, re, t ); // re = |z|
  mad_tpsa_atan2(im, re, im); // im = arg(z)
  mad_ctpsa_cplx(t , im, c );

  REL_TMPR(re); REL_TMPR(im); REL_TMPR(t);

  DBGTPSA(c); DBGFUN(<-);
}

log_t mad_ctpsa_equt (const ctpsa_t *a, const tpsa_t *b, num_t tol)
{
  assert(a && b); DBGFUN(->);
  ensure(a->d == b->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  log_t res = mad_ctpsa_equ(a, t, tol);
  REL_TMPC(t);
  DBGFUN(<-); return res;
}

void mad_ctpsa_addt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_add (a, t, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_subt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_sub (a, t, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_tsub (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_sub (t, b, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_mult (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_mul (a, t, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_divt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_div (a, t, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_tdiv (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_div (t, b, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_powt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_pow (a, t, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_tpow (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_pow (t, b, c);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_poisst(const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c, int nv)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx   (b, NULL, t);
  mad_ctpsa_poisson(a, t, c, nv);
  REL_TMPC(t);
  DBGFUN(<-);
}

void mad_ctpsa_tpoiss(const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx   (a, NULL, t);
  mad_ctpsa_poisson(t, b, c, nv);
  REL_TMPC(t);
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
