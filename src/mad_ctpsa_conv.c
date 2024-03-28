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

#include "mad_mem.h"
#include "mad_ctpsa_impl.h"

// --- conversion (special cases) ---------------------------------------------o

void
mad_ctpsa_real (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c); DBGFUN(->);
  if (DEBUG > 1) mad_ctpsa_debug(t,"t",__func__,__LINE__,0);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  FUN(copy0)((const tpsa_t*)t, c);

  c->coef[0] = creal(t->coef[0]);

  if (!c->nz) { FUN(setval)(c, c->coef[0]); DBGFUN(<-); return; }

  bit_t cnz = c->nz;
  TPSA_SCAN_Z(c) {
    log_t nz = 0;
    TPSA_SCAN_O(c) c->coef[i] = creal(t->coef[i]), nz |= !!c->coef[i];
    if (!nz) cnz = mad_bit_clr(cnz,o);
  }
  if (c->nz != cnz) c->nz = cnz, FUN(adjust0)(c);
  DBGTPSA(c); DBGFUN(<-);
}

void
mad_ctpsa_imag (const ctpsa_t *t, tpsa_t *c)
{
  assert(t && c); DBGFUN(->);
  if (DEBUG > 1) mad_ctpsa_debug(t,"t",__func__,__LINE__,0);
  ensure(t->d == c->d, "incompatibles GTPSA (descriptors differ)");

  FUN(copy0)((const tpsa_t*)t, c);

  c->coef[0] = cimag(t->coef[0]);

  if (!c->nz) { FUN(setval)(c, c->coef[0]); DBGFUN(<-); return; }

  bit_t cnz = c->nz;
  TPSA_SCAN_Z(c) {
    log_t nz = 0;
    TPSA_SCAN_O(c) c->coef[i] = cimag(t->coef[i]), nz |= !!c->coef[i];
    if (!nz) cnz = mad_bit_clr(cnz,o);
  }
  if (c->nz != cnz) c->nz = cnz, FUN(adjust0)(c);
  DBGTPSA(c); DBGFUN(<-);
}

void
mad_ctpsa_cplx (const tpsa_t *re_, const tpsa_t *im_, ctpsa_t *c)
{
  assert((re_ || im_) && c); DBGFUN(->);
  const tpsa_t *re = re_ ? re_ : im_; DBGTPSA(re);
  const tpsa_t *im = im_ ? im_ : re_; DBGTPSA(im);
  ensure(re->d == c->d && im->d == c->d, "incompatibles GTPSA (descriptors differ)");

  FUN(copy00)(re, im, (tpsa_t*)c);

  c->coef[0] = (re_ ? re_->coef[0] : 0) + (im_ ? im_->coef[0] : 0)*I;

  if (!c->nz) { mad_ctpsa_setval(c, c->coef[0]); DBGFUN(<-); return; }

  switch(!!re_ + 2*!!im_) {
  case 1: { TPSA_SCAN(c) c->coef[i] = re_->coef[i];   break; }
  case 2: { TPSA_SCAN(c) c->coef[i] = im_->coef[i]*I; break; }
  case 3: {
    TPSA_SCAN(c) {
      c->coef[i] = 0;
      if (mad_bit_tst(re_->nz,o)) c->coef[i] += re_->coef[i];
      if (mad_bit_tst(im_->nz,o)) c->coef[i] += im_->coef[i]*I;
    }}
  }
  if (DEBUG > 1) mad_ctpsa_debug(c,"c",__func__,__LINE__,0);
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
  REL_TMPR(im); REL_TMPR(re); DBGFUN(<-);
}

void mad_ctpsa_carg (const ctpsa_t *a, tpsa_t *c)
{
  assert(a && c); DBGFUN(->);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  tpsa_t *re = GET_TMPR(a), *im = GET_TMPR(a);
  mad_ctpsa_real(a, re);
  mad_ctpsa_imag(a, im);
  mad_tpsa_atan2(im, re, c);
  REL_TMPR(im); REL_TMPR(re); DBGFUN(<-);
}

void mad_ctpsa_rect (const ctpsa_t *a, ctpsa_t *c)
{
  assert(a && c); DBGFUN(->);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  tpsa_t *re = GET_TMPR(c), *im = GET_TMPR(c),
         *st = GET_TMPR(c), *ct = GET_TMPR(c);
  mad_ctpsa_real (a , re); // rho
  mad_ctpsa_imag (a , im); // arg
  mad_tpsa_sincos(im, st, ct);
  mad_tpsa_mul   (re, st, im); // im = rho * sin(arg)
  mad_tpsa_mul   (re, ct, st); // re = rho * cos(arg)
  mad_ctpsa_cplx (st, im, c );
  REL_TMPR(ct); REL_TMPR(st); REL_TMPR(im); REL_TMPR(re); DBGFUN(<-);
}

void mad_ctpsa_polar (const ctpsa_t *a, ctpsa_t *c)
{
  assert(a && c); DBGFUN(->);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  tpsa_t *re = GET_TMPR(c), *im = GET_TMPR(c), *t = GET_TMPR(c);
  mad_ctpsa_real(a , re);
  mad_ctpsa_imag(a , im);
  mad_tpsa_hypot(im, re, t ); // re = |z|
  mad_tpsa_atan2(im, re, im); // im = arg(z)
  mad_ctpsa_cplx(t , im, c );
  REL_TMPR(t); REL_TMPR(im); REL_TMPR(re); DBGFUN(<-);
}

log_t mad_ctpsa_equt (const ctpsa_t *a, const tpsa_t *b, num_t tol_)
{
  assert(a && b); DBGFUN(->);
  ensure(a->d == b->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  log_t res = mad_ctpsa_equ(a, t, tol_);
  REL_TMPC(t); DBGFUN(<-); return res;
}

void mad_ctpsa_dift (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_dif (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_tdif (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_dif (t, b, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_addt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_add (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_subt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_sub (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_tsub (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_sub (t, b, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_mult (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_mul (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_divt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_div (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_tdiv (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_div (t, b, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_powt (const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx(b, NULL, t);
  mad_ctpsa_pow (a, t, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_tpow (const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx(a, NULL, t);
  mad_ctpsa_pow (t, b, c);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_poisbrat(const ctpsa_t *a, const tpsa_t *b, ctpsa_t *c, int nv)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(b);
  mad_ctpsa_cplx   (b, NULL, t);
  mad_ctpsa_poisbra(a, t, c, nv);
  REL_TMPC(t); DBGFUN(<-);
}

void mad_ctpsa_tpoisbra(const tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv)
{
  assert(a && b && c); DBGFUN(->);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ctpsa_t *t = GET_TMPC(a);
  mad_ctpsa_cplx   (a, NULL, t);
  mad_ctpsa_poisbra(t, b, c, nv);
  REL_TMPC(t); DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
