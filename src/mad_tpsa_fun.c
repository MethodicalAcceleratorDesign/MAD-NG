/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA functions module implementation
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

#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

#include "mad_cst.h"
#include "mad_log.h"
#include "mad_num.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#include <complex.h>
#endif

// --- local ------------------------------------------------------------------o

enum { MANUAL_EXPANSION_ORD = 5 };

static inline void
fixed_point_iteration(const T *a, T *c, ord_t iter, NUM ord_coef[iter+1])
{
  assert(a && c && ord_coef);
  assert(iter >= 1); // ord 0 treated outside

  T *acp = GET_TMPX(a);
  if (iter >=2) FUN(copy)(a,acp); // save copy before scaling for aliasing

  // iter 1
  FUN(scl)(a, ord_coef[1], c);
  FUN(set0)(c, 0, ord_coef[0]);

  // iter 2..iter
  if (iter >= 2) {
    T *pow = GET_TMPX(a);
    T *tmp = GET_TMPX(a), *t;
    FUN(set0)(acp,0,0);  // clear a0
    FUN(copy)(acp,pow);  // already did ord 1

    for (ord_t i = 2; i <= iter; ++i) {
      FUN(mul)(acp,pow,tmp); // never use t[0]!
      FUN(acc)(tmp,ord_coef[i],c);
      SWAP(pow,tmp,t);
    }
    REL_TMPX(tmp), REL_TMPX(pow);
  }
  REL_TMPX(acp);
}

static inline void
sincos_fixed_point(const T *a, T *s, T *c, ord_t iter_s,
                   NUM sin_coef[iter_s+1], ord_t iter_c, NUM cos_coef[iter_c+1])
{
  assert(a && s && c && sin_coef && cos_coef);
  assert(iter_s >= 1 && iter_c >= 1);  // ord 0 treated outside

  T *acp = GET_TMPX(a);
  ord_t iter = MAX(iter_s,iter_c);
  if (iter >= 2) FUN(copy)(a,acp); // save copy before scaling for aliasing

  // iter 1
  FUN(scl)(a, sin_coef[1], s); FUN(set0)(s, 0, sin_coef[0]);
  FUN(scl)(a, cos_coef[1], c); FUN(set0)(c, 0, cos_coef[0]);

  if (iter >= 2) {
    T *pow = GET_TMPX(a);
    T *tmp = GET_TMPX(a), *t;
    FUN(set0)(acp,0,0);
    FUN(copy)(acp,pow);

    for (ord_t i = 2; i <= iter; ++i) {
      FUN(mul)(acp,pow,tmp); // never use t[0]!
      if (i <= iter_s) FUN(acc)(tmp,sin_coef[i],s);
      if (i <= iter_c) FUN(acc)(tmp,cos_coef[i],c);
      SWAP(pow,tmp,t);
    }
    REL_TMPX(tmp), REL_TMPX(pow);
  }
  REL_TMPX(acp);
}

// --- public -----------------------------------------------------------------o

void
FUN(inv) (const T *a, NUM v, T *c) // c = v/a
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ensure(a->coef[0] != 0, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, v/a->coef[0]); return; }

  NUM ord_coef[to+1], a0 = a->coef[0], _a0 = 1 / a0;
  ord_coef[0] = _a0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0;

  fixed_point_iteration(a,c,to,ord_coef);
  FUN(scl)(c,v,c);
}

void
FUN(sqrt) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = sqrt(a->coef[0]);
  if (errno) error(strerror(errno));
  // SELECT(ensure(a->coef[0] >= 0, "invalid domain"), );

  if (a->coef[0] == 0) { FUN(clear)(c); return; }

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1], a0 = a->coef[0], _a0 = 1 / a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-3);

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(invsqrt) (const T *a, NUM v, T *c)  // v/sqrt(a)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = 1/sqrt(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(a->coef[0] SELECT(> 0, != 0), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, v*f0); return; }

  NUM ord_coef[to+1], a0 = a->coef[0], _a0 = 1 / a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-1);

  fixed_point_iteration(a,c,to,ord_coef);
  FUN(scl)(c,v,c);
}

void
FUN(exp) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = exp(a->coef[0]);
  // if (errno) error(strerror(errno));

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1]; //, a0 = a->coef[0];
  ord_coef[0] = f0;
  for (int o = 1; o <= to; ++o)
    ord_coef[o] = ord_coef[o-1] / o;

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(log) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = log(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(a->coef[0] SELECT(> 0, != 0), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1], a0 = a->coef[0], _a0 = 1 / a0;
  ord_coef[0] = f0;
  ord_coef[1] = _a0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / o * (o-1);

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(sincos) (const T *a, T *s, T *c)
{
  assert(a && s && c);
  ensure(a->d == s->d && a->d == c->d, "incompatible GTPSA (descriptors differ)");

  ord_t sto = MIN(s->mo,s->d->trunc),
        cto = MIN(c->mo,c->d->trunc);

  NUM a0 = a->coef[0], sa = sin(a0), ca = cos(a0);

  if (a->hi == 0) {
    FUN(scalar)(s, sa);
    FUN(scalar)(c, ca);
    return;
  }

  if (!sto || !cto) {
    if (!sto) FUN(scalar)(s, sa);
    else      FUN(sin)(a,s);
    if (!cto) FUN(scalar)(c, ca);
    else      FUN(cos)(a,c);
    return;
  }

  // ord 0, 1
  NUM sin_coef[sto+1], cos_coef[cto+1];
  sin_coef[0] = sa;  cos_coef[0] =  ca;
  sin_coef[1] = ca;  cos_coef[1] = -sa;

  // ords 2..to
  for (ord_t o = 2; o <= sto; ++o )
    sin_coef[o] = -sin_coef[o-2] / (o*(o-1));
  for (ord_t o = 2; o <= cto; ++o )
    cos_coef[o] = -cos_coef[o-2] / (o*(o-1));

  sincos_fixed_point(a,s,c, sto,sin_coef, cto,cos_coef);
}

void
FUN(sin) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, sin(a->coef[0])); return; }

  NUM ord_coef[to+1], a0 = a->coef[0];
  ord_coef[0] = sin(a0);
  ord_coef[1] = cos(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(cos) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, cos(a->coef[0])); return; }

  NUM ord_coef[to+1], a0 = a->coef[0];
  ord_coef[0] =  cos(a0);
  ord_coef[1] = -sin(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(tan) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = tan(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(ca != 0, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM f2 = f0*f0;
  switch(to) {
  case 5: ord_coef[5] = 2.0/15 + f2*(17.0/15 + f2*(2 + f2)); /* FALLTHRU */
  case 4: ord_coef[4] = f0*(2.0/3 + f2*(5.0/3 + f2));        /* FALLTHRU */
  case 3: ord_coef[3] = 1.0/3 + f2*(4.0/3 + f2);             /* FALLTHRU */
  case 2: ord_coef[2] = f0*(1 + f2);                         /* FALLTHRU */
  case 1: ord_coef[1] = 1 + f2;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(cot) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = tan(M_PI_2 - a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(sa != 0, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], sa = sin(a0), ca = cos(a0);
  NUM xc1 = 1/sa, xc2 = xc1*xc1, sa2 = sa*sa, ca2 = ca*ca;
  switch(to) {
  case 5: ord_coef[5] = -(2*sa2 + 3*sa2*ca2 + 10*ca2 + 5*ca2*ca2) *xc2*xc2*xc2 /15; /* FALLTHRU */
  case 4: ord_coef[4] =  (2*ca  +   ca2*ca) *xc2*xc2*xc1 /3; /* FALLTHRU */
  case 3: ord_coef[3] = -(  sa2 + 3*ca2   ) *xc2*xc2 /3;     /* FALLTHRU */
  case 2: ord_coef[2] = ca                  *xc2*xc1;        /* FALLTHRU */
  case 1: ord_coef[1] = -1                  *xc2;            /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(sinc) (const T *a, T *c)
{
// SIN(P)/P = 1 - P^2/3! + P^4/5! - P^6/7! + ...
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

#ifdef MAD_CTPSA_IMPL
  NUM f0 = mad_num_sinc (a->coef[0]);
#else
  NUM f0 = mad_cnum_sinc(a->coef[0]);
#endif

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0; // was 1, certainly buggy, see sin, TODO check!!!
  ord_coef[1] = 0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o * (o+1));

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(sincosh) (const T *a, T *sh, T *ch)
{
  assert(a && sh && ch);
  ensure(a->d == sh->d && a->d == ch->d, "incompatible GTPSA (descriptors differ)");

  ord_t sto = MIN(sh->mo,sh->d->trunc),
        cto = MIN(ch->mo,ch->d->trunc);

  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);
  if (a->hi == 0) {
    FUN(scalar)(sh, sa);
    FUN(scalar)(ch, ca);
    return;
  }

  if (!sto || !cto) {
    if (!sto) FUN(scalar)(sh, sa);
    else      FUN(sinh)(a,sh);
    if (!cto) FUN(scalar)(ch, ca);
    else      FUN(cosh)(a,ch);
    return;
  }

  // ord 0, 1
  NUM sin_coef[sto+1], cos_coef[cto+1];
  sin_coef[0] = sa;  cos_coef[0] = ca;
  sin_coef[1] = ca;  cos_coef[1] = sa;

  // ords 2..to
  for (int o = 2; o <= sto; ++o )
    sin_coef[o] = sin_coef[o-2] / (o*(o-1));
  for (int o = 2; o <= cto; ++o )
    cos_coef[o] = cos_coef[o-2] / (o*(o-1));

  sincos_fixed_point(a,sh,ch, sto,sin_coef, cto,cos_coef);
}

void
FUN(sinh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = sinh(a->coef[0]);
  if (errno) error(strerror(errno));

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1], a0 = a->coef[0];
  ord_coef[0] = f0;
  ord_coef[1] = cosh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(cosh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = cosh(a->coef[0]);
  if (errno) error(strerror(errno));

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  NUM ord_coef[to+1], a0 = a->coef[0];
  ord_coef[0] = f0;
  ord_coef[1] = sinh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(tanh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = tanh(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(ca != 0, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);
  NUM xc1 = 1/ca, xc2 = xc1*xc1, sa2 = sa*sa, ca2 = ca*ca;
  switch(to) {
  case 5: ord_coef[5] = (2*ca2 - 3*ca2*sa2 - 10*sa2 + 5*sa2*sa2) *xc2*xc2*xc2 /15; /* FALLTHRU */
  case 4: ord_coef[4] = (2*sa  -   sa2*sa) *xc2*xc2*xc1 /3; /* FALLTHRU */
  case 3: ord_coef[3] = (- ca2 + 3*sa2   ) *xc2*xc2 /3;     /* FALLTHRU */
  case 2: ord_coef[2] = -sa                *xc2*xc1;        /* FALLTHRU */
  case 1: ord_coef[1] =  1                 *xc2;            /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(coth) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = tanh(a->coef[0]);
  if (errno) error(strerror(errno));
  ensure(f0 != 0, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, 1/f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);
  NUM xc1 = 1/sa, xc2 = xc1*xc1, sa2 = sa*sa, ca2 = ca*ca;
  switch(to) {
  case 5: ord_coef[5] = (2*sa2 + 3*sa2*ca2 - 10*ca2 - 5*ca2*ca2) *xc2*xc2*xc2 /15; /* FALLTHRU */
  case 4: ord_coef[4] = (2*ca  +   ca2*ca) *xc2*xc2*xc1 /3; /* FALLTHRU */
  case 3: ord_coef[3] = (  sa2 - 3*ca2   ) *xc2*xc2 /3;     /* FALLTHRU */
  case 2: ord_coef[2] = ca                 *xc2*xc1;        /* FALLTHRU */
  case 1: ord_coef[1] = -1                 *xc2;            /* FALLTHRU */
  case 0: ord_coef[0] = ca                 *xc1;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(asin) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = asin(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(SELECT(fabs(a0) < 1, a2 != 1), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // asin(x) = -i*ln(i*x + sqrt(1-x^2))
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_scl(c, -I, c);
#else
    ctpsa_t *t = mad_ctpsa_newd(c->d, to);
    mad_ctpsa_complex(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_scl(t, -I, t);
    mad_ctpsa_real(t, c);
    mad_ctpsa_del(t);
#endif
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / sqrt(1 - a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = (3    + 24*a2 + 8*a2*a2) *xc2*xc2*xc2*xc2*xc1 /40; /* FALLTHRU */
  case 4: ord_coef[4] = (3*a0 +  2*a2*a0) *xc2*xc2*xc2*xc1 /8; /* FALLTHRU */
  case 3: ord_coef[3] = (1    +  2*a2   ) *xc2*xc2*xc1 /6;     /* FALLTHRU */
  case 2: ord_coef[2] =               a0  *xc2*xc1 /2;         /* FALLTHRU */
  case 1: ord_coef[1] =                    xc1;                /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(acos) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = acos(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(SELECT(fabs(a0) < 1, a2 != 1), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acos(x) = -i*ln(x+i*sqrt(1-x^2)) = -asin(x)+pi/2
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_axpb(I, c, M_PI_2, c);
#else
    ctpsa_t *t = mad_ctpsa_newd(c->d, to);
    mad_ctpsa_complex(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_axpb(I, t, M_PI_2, t);
    mad_ctpsa_real(t, c);
    mad_ctpsa_del(t);
#endif
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / sqrt(1 - a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = -(3    + 24*a2 + 8*a2*a2) *xc2*xc2*xc2*xc2*xc1 /40; /* FALLTHRU */
  case 4: ord_coef[4] = -(3*a0 +  2*a2*a0) *xc2*xc2*xc2*xc1 /8; /* FALLTHRU */
  case 3: ord_coef[3] = -(1    +  2*a2   ) *xc2*xc2*xc1 /6;     /* FALLTHRU */
  case 2: ord_coef[2] = -              a0  *xc2*xc1 /2;         /* FALLTHRU */
  case 1: ord_coef[1] = -                   xc1;                /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(atan) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = atan(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(a2 != -1, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // atan(x) = i/2 ln((i+x) / (i-x))
    ctpsa_t *tn = mad_ctpsa_newd(c->d, to);
    ctpsa_t *td = mad_ctpsa_newd(c->d, to);
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_copy(a, tn);
    mad_ctpsa_set0(tn, 1, I);
    mad_ctpsa_copy(a, td);
    mad_ctpsa_set0(td, -1, I);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, I/2, c);
#else
    ctpsa_t *t = mad_ctpsa_newd(c->d, to);
    mad_ctpsa_complex(a, NULL, tn);
    mad_ctpsa_set0(tn, 1, I);
    mad_ctpsa_complex(a, NULL, td);
    mad_ctpsa_set0(td, -1, I);
    mad_ctpsa_logxdy(tn, td, t);
    mad_ctpsa_scl(t, I/2, t);
    mad_ctpsa_real(t, c);
    mad_ctpsa_del(t);
#endif
    mad_ctpsa_del(td);
    mad_ctpsa_del(tn);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / (1 + a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] =  (1.0/5 + a2*a2 - 2*a2) * xc2*xc2*xc1; /* FALLTHRU */
  case 4: ord_coef[4] =  (a0    - a2*a0)        * xc2*xc2;     /* FALLTHRU */
  case 3: ord_coef[3] = -(1.0/3 - a2)           * xc2*xc1;     /* FALLTHRU */
  case 2: ord_coef[2] = -a0                     * xc2;         /* FALLTHRU */
  case 1: ord_coef[1] =                           xc1;         /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(acot) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = atan(1/a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(a2 != -1, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acot(x) = i/2 ln((x-i) / (x+i))
    ctpsa_t *tn = mad_ctpsa_newd(c->d, to);
    ctpsa_t *td = mad_ctpsa_newd(c->d, to);
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_copy(a, tn);
    mad_ctpsa_set0(tn, -I, 1);
    mad_ctpsa_copy(a, td);
    mad_ctpsa_set0(td, I, 1);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, I/2, c);
#else
    ctpsa_t *t = mad_ctpsa_newd(c->d, to);
    mad_ctpsa_complex(a, NULL, tn);
    mad_ctpsa_set0(tn, -I, 1);
    mad_ctpsa_complex(a, NULL, td);
    mad_ctpsa_set0(td, I, 1);
    mad_ctpsa_logxdy(tn, td, t);
    mad_ctpsa_scl(t, I/2, t);
    mad_ctpsa_real(t, c);
    mad_ctpsa_del(t);
#endif
    mad_ctpsa_del(td);
    mad_ctpsa_del(tn);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / (1 + a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = -(1.0/5 + a2*a2 - 2*a2) * xc2*xc2*xc1; /* FALLTHRU */
  case 4: ord_coef[4] = -(a0    - a2*a0)        * xc2*xc2;     /* FALLTHRU */
  case 3: ord_coef[3] =  (1.0/3 - a2)           * xc2*xc1;     /* FALLTHRU */
  case 2: ord_coef[2] =           a0            * xc2;         /* FALLTHRU */
  case 1: ord_coef[1] = -                         xc1;         /* FALLTHRU */
  case 0: ord_coef[0] = f0; // was M_PI_2 - atan(a0);
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(asinh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = asinh(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(a2 != -1, "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // asinh(x) = log(x + sqrt(x^2+1))
    FUN(logaxpsqrtbpcx2)(a, 1, 1, 1, c);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / sqrt(1 + a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = (3    - 24*a2 + 8*a2*a2) *xc2*xc2*xc2*xc2*xc1 /40; /* FALLTHRU */
  case 4: ord_coef[4] = (3*a0 -  2*a2*a0) * xc2*xc2*xc2*xc1 /8; /* FALLTHRU */
  case 3: ord_coef[3] = (-1   +  2*a2   ) * xc2*xc2*xc1 /6;     /* FALLTHRU */
  case 2: ord_coef[2] = -          a0     * xc2*xc1 /2;         /* FALLTHRU */
  case 1: ord_coef[1] =                     xc1;                /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(acosh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = acosh(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(SELECT(a0 > 1, a2 != 1), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acosh(x) = ln(x + sqrt(x^2-1))
    FUN(logaxpsqrtbpcx2)(a, 1, -1, 1, c);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / sqrt(a2 - 1), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = (3     + 24*a2 + 8*a2*a2) *xc2*xc2*xc2*xc2*xc1 /40; /* FALLTHRU */
  case 4: ord_coef[4] = (-3*a0 -  2*a2*a0) *xc2*xc2*xc2*xc1 /8; /* FALLTHRU */
  case 3: ord_coef[3] = (1     +  2*a2   ) *xc2*xc2*xc1 /6;     /* FALLTHRU */
  case 2: ord_coef[2] = -           a0     *xc2*xc1 /2;         /* FALLTHRU */
  case 1: ord_coef[1] =                     xc1;                /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(atanh) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = atanh(a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(SELECT(fabs(a0) < 1, a2 != 1), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // atanh(x) = 1/2 ln((1+x) / (1-x))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(set0)(tn, 1, 1);
    FUN(copy)(a, td);
    FUN(set0)(td, 1, -1);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / (1 - a2), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = (1.0/5 + a2*a2 + 2*a2) * xc2*xc2*xc1; /* FALLTHRU */
  case 4: ord_coef[4] = (a0    + a2*a0)        * xc2*xc2;     /* FALLTHRU */
  case 3: ord_coef[3] = (1.0/3 + a2)           * xc2*xc1;     /* FALLTHRU */
  case 2: ord_coef[2] =  a0                    * xc2;         /* FALLTHRU */
  case 1: ord_coef[1] =                          xc1;         /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(acoth) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM f0 = atanh(1/a->coef[0]);
  if (errno) error(strerror(errno));
  // ensure(SELECT(fabs(a0) > 1, a0 != 0 && a2 != 1), "invalid domain");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (!to || a->hi == 0) { FUN(scalar)(c, f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acoth(x) = 1/2 ln((x+1) / (x-1))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(set0)(tn, 1, 1);
    FUN(copy)(a, td);
    FUN(set0)(td, -1, 1);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn);
    return;
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;
  NUM xc1 = 1 / (a2 - 1), xc2 = xc1*xc1;
  switch(to) {
  case 5: ord_coef[5] = (-1.0/5 - a2*a2 - 2*a2) *xc2*xc2*xc1; /* FALLTHRU */
  case 4: ord_coef[4] = (a0     + a2*a0)        *xc2*xc2;     /* FALLTHRU */
  case 3: ord_coef[3] = (-1.0/3 - a2   )        *xc2*xc1;     /* FALLTHRU */
  case 2: ord_coef[2] =              a0         *xc2;         /* FALLTHRU */
  case 1: ord_coef[1] = -                        xc1;         /* FALLTHRU */
  case 0: ord_coef[0] = f0;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

void
FUN(erf) (const T *a, T *c)
{
  // ERF(X) is the integral from 0 to x from [2/sqrt(pi) * exp(-x*x)]
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  ord_t to = MIN(c->mo,c->d->trunc);
  if (to > MANUAL_EXPANSION_ORD) {
    ensure(to <= MANUAL_EXPANSION_ORD, "NYI");
    // TODO: use Faddeva function from MIT??
  }

  NUM ord_coef[MANUAL_EXPANSION_ORD+1];
  NUM a0 = a->coef[0], a2 = a0*a0;

  static const // coefs from Mathematica
  num_t b1 =  0.254829592,
        b2 = -0.284496736,
        b3 =  1.421413741,
        b4 = -1.453152027,
        b5 =  1.061405429,
        p  =  0.327591100,
        p2 =  0.886226925452758013649083741671; // sqrt(atan(1)),
  NUM   t  = 1 / (1 + p*a0),
        e1 = exp(-a2),
        e2 = 1 - t*(b1+t*(b2+t*(b3+t*(b4+t*b5))))*e1;

  switch(to) {
  case 5: ord_coef[5] = (16*a2*a2 - 48*a2 + 12) / 120*e1 / p2; /* FALLTHRU */
  case 4: ord_coef[4] = (12*a0    -  8*a2*a0)   /  24*e1 / p2; /* FALLTHRU */
  case 3: ord_coef[3] = (-1       +  2*a2)      /   3*e1 / p2; /* FALLTHRU */
  case 2: ord_coef[2] = -                          a0*e1 / p2; /* FALLTHRU */
  case 1: ord_coef[1] =                               e1 / p2; /* FALLTHRU */
  case 0: ord_coef[0] =                               e2;
  }

  fixed_point_iteration(a,c,to,ord_coef);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(inv_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(inv)(a, CNUM(v), c); }

void FUN(invsqrt_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(invsqrt)(a, CNUM(v), c); }

#endif

// --- end --------------------------------------------------------------------o
