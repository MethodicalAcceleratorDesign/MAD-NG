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
#define   MAD_TPSA_NOHELPER
#include "mad_tpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#define   MAD_TPSA_NOHELPER
#include "mad_ctpsa_impl.h"
#endif

// --- local ------------------------------------------------------------------o

enum { MANUAL_EXPANSION_ORD = 6 };

static inline void
fun_fix_point (const T *a, T *c, ord_t iter, NUM ord_coef[iter+1])
{
  assert(a && c && ord_coef);
  assert(iter >= 1); // ord 0 treated outside

  T *acp = GET_TMPX(c);
  if (iter >=2) FUN(copy)(a,acp); // make copy before scaling for aliasing

  // iter 1
  FUN(scl)(a, ord_coef[1], c);
  FUN(set0)(c, 0, ord_coef[0]);

  // iter 2..iter
  if (iter >= 2) {
    T *pow = GET_TMPX(c);
    T *tmp = GET_TMPX(c), *t;
    FUN(set0)(acp,0,0);  // clear a0
    FUN(copy)(acp,pow);  // already did ord 1

    for (ord_t i = 2; i <= iter; ++i) {
      FUN(mul)(acp,pow,tmp);
      FUN(acc)(tmp,ord_coef[i],c);
      SWAP(pow,tmp,t);
    }

    if ((iter-1) & 1) SWAP(pow,tmp,t); // enforce even number of swaps
    REL_TMPX(tmp), REL_TMPX(pow);
  }
  REL_TMPX(acp);
}

static inline void
sincos_fix_point (const T *a, T *s, T *c, ord_t iter_s,
                  NUM sin_coef[iter_s+1], ord_t iter_c, NUM cos_coef[iter_c+1])
{
  assert(a && s && c && sin_coef && cos_coef);
  assert(iter_s >= 1 && iter_c >= 1);  // ord 0 treated outside

  T *acp = GET_TMPX(c);
  ord_t iter = MAX(iter_s,iter_c);
  if (iter >= 2) FUN(copy)(a,acp); // make copy before scaling for aliasing

  // iter 1
  FUN(scl)(a, sin_coef[1], s); FUN(set0)(s, 0, sin_coef[0]);
  FUN(scl)(a, cos_coef[1], c); FUN(set0)(c, 0, cos_coef[0]);

  if (iter >= 2) {
    T *pow = GET_TMPX(c);
    T *tmp = GET_TMPX(c), *t;
    FUN(set0)(acp,0,0);
    FUN(copy)(acp,pow);

    for (ord_t i = 2; i <= iter; ++i) {
      FUN(mul)(acp,pow,tmp);
      if (i <= iter_s) FUN(acc)(tmp,sin_coef[i],s);
      if (i <= iter_c) FUN(acc)(tmp,cos_coef[i],c);
      SWAP(pow,tmp,t);
    }

    if ((iter-1) & 1) SWAP(pow,tmp,t); // enforce even number of swaps
    REL_TMPX(tmp), REL_TMPX(pow);
  }
  REL_TMPX(acp);
}

// --- public -----------------------------------------------------------------o

void
FUN(inv) (const T *a, NUM v, T *c) // c = v/a    // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain inv("FMT")", VAL(a0));
  NUM f0 = 1/a0;

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,v*f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * f0;

  fun_fix_point(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
}

void
FUN(invsqrt) (const T *a, NUM v, T *c) // v/sqrt(a),checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain invsqrt("FMT")", VAL(a0));
  NUM f0 = 1/sqrt(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,v*f0); return; }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-1);

  fun_fix_point(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
}

void
FUN(sqrt) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain sqrt("FMT")", VAL(a0));
  NUM f0 = sqrt(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-3);

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(exp) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = exp(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1]; //, a0 = a->coef[0];
  ord_coef[0] = f0;
  for (int o = 1; o <= to; ++o)
    ord_coef[o] = ord_coef[o-1] / o;

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(log) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain log("FMT")", VAL(a0));
  NUM f0 = log(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  ord_coef[1] = _a0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / o * (o-1);

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(pow) (const T *a, const T *b, T *c)          // checked for real and complex
{
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(mul)(b,t,c);
  FUN(exp)(c,c);
  REL_TMPX(t);
}

void
FUN(pown) (const T *a, NUM v, T *c)              // checked for real and complex
{
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(scl)(t,v,c);
  FUN(exp)(c,c);
  REL_TMPX(t);
}

void
FUN(sincos) (const T *a, T *s, T *c)             // checked for real and complex
{
  assert(a && s && c);
  ensure(a->d == s->d && a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sin(a0), ca = cos(a0);

  if (a->hi == 0) {
    FUN(scalar)(s, sa);
    FUN(scalar)(c, ca);
    return;
  }

  ord_t sto = MIN(s->mo,s->d->to),
        cto = MIN(c->mo,c->d->to);
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

  sincos_fix_point(a,s,c, sto,sin_coef, cto,cos_coef);
}

void
FUN(sin) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sin(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cos(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(cos) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cos(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = -sin(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(tan) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(cos(a0) != 0, "invalid domain tan("FMT")", VAL(a0));
  NUM f0 = tan(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[to+1], f2 = f0*f0;
  switch(to) {
  case 6: ord_coef[6] = f0*(17./45 + f2*(77./45 + f2*(7./3 + f2))); /* FALLTHRU */
  case 5: ord_coef[5] = 2./15 + f2*(17./15 + f2*(2 + f2));          /* FALLTHRU */
  case 4: ord_coef[4] = f0*(2./3 + f2*(5./3 + f2));                 /* FALLTHRU */
  case 3: ord_coef[3] = 1./3 + f2*(4./3 + f2);                      /* FALLTHRU */
  case 2: ord_coef[2] = f0*(1 + f2);                                /* FALLTHRU */
  case 1: ord_coef[1] = 1 + f2;                                     /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                         break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(cot) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(sin(a0) != 0, "invalid domain cot("FMT")", VAL(a0));
  NUM f0 = tan(M_PI_2 - a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[to+1], f2 = f0*f0;
  switch(to) {
  case 6: ord_coef[6] = f0*(17./45 + f2*(77./45 + f2*(7./3 + f2))); /* FALLTHRU */
  case 5: ord_coef[5] = -2./15 - f2*(17./15 + f2*(2 + f2));         /* FALLTHRU */
  case 4: ord_coef[4] = f0*(2./3 + f2*(5./3 + f2));                 /* FALLTHRU */
  case 3: ord_coef[3] = -1./3 - f2*(4./3 + f2);                     /* FALLTHRU */
  case 2: ord_coef[2] = f0*(1 + f2);                                /* FALLTHRU */
  case 1: ord_coef[1] = -1 - f2;                                    /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                         break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(sinc) (const T *a, T *c)                 // NOT checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  // With IEEE 754: cos(1e-10) = 1, sin(1e-10) = 0
  NUM a0 = a->coef[0];
#ifdef MAD_CTPSA_IMPL
  cnum_t f0 = mad_cnum_sinc(a0);
#else
  num_t  f0 = mad_num_sinc (a0);
#endif

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1];

  if (fabs(f0) < 1 - 1e-10) {
    T *t = GET_TMPX(c);
    FUN(sin)(a,t);
    FUN(div)(t,a,c);
    REL_TMPX(t);
    return;
  } else
  // prefer explicit above
  if (fabs(f0) < 1 - 1e-10) {
    NUM sa = -sin(a0), ca = cos(a0), _a0 = 1/a0, term;
    num_t fo = 1;
    ord_coef[0] = f0;
    for (int o = 1; o <= to; ++o) {
      fo *= o;
      term = o & 1 ? (ca=-ca,ca) : (sa=-sa,sa);
      ord_coef[o] = -(ord_coef[o-1] + term/fo)*_a0;
    }
  } else {  // |a0| < 1e-10
    ord_coef[0] = 1;
    ord_coef[1] = 0;
    for (int o = 2; o <= to; ++o)
      ord_coef[o] = -ord_coef[o-2] / (o * (o+1));
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(sincosh) (const T *a, T *sh, T *ch)          // checked for real and complex
{
  assert(a && sh && ch);
  ensure(a->d == sh->d && a->d == ch->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);

  if (a->hi == 0) {
    FUN(scalar)(sh, sa);
    FUN(scalar)(ch, ca);
    return;
  }

  ord_t sto = MIN(sh->mo,sh->d->to),
        cto = MIN(ch->mo,ch->d->to);
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

  sincos_fix_point(a,sh,ch, sto,sin_coef, cto,cos_coef);
}

void
FUN(sinh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sinh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cosh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(cosh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cosh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = sinh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(tanh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[to+1], f2 = f0*f0;
  switch(to) {
  case 6: ord_coef[6] = f0*(-17./45 + f2*(77./45 + f2*(-7./3 + f2))); /* FALLTHRU */
  case 5: ord_coef[5] = 2./15 + f2*(-17./15 + f2*(2 - f2));           /* FALLTHRU */
  case 4: ord_coef[4] = f0*(2./3 + f2*(-5./3 + f2));                  /* FALLTHRU */
  case 3: ord_coef[3] = -1./3 + f2*(4./3 - f2);                       /* FALLTHRU */
  case 2: ord_coef[2] = f0*(-1 + f2);                                 /* FALLTHRU */
  case 1: ord_coef[1] = 1 - f2;                                       /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                           break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(coth) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);
  ensure(f0 != 0, "invalid domain coth("FMT")", VAL(a0));
  f0 = 1/f0;

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t);
    return;
  }

  NUM ord_coef[to+1], f2 = f0*f0;
  switch(to) {
  case 6: ord_coef[6] = f0*(-17./45 + f2*(77./45 + f2*(-7./3 + f2))); /* FALLTHRU */
  case 5: ord_coef[5] = 2./15 + f2*(-17./15 + f2*(2 - f2));           /* FALLTHRU */
  case 4: ord_coef[4] = f0*(2./3 + f2*(-5./3 + f2));                  /* FALLTHRU */
  case 3: ord_coef[3] = -1./3 + f2*(4./3 - f2);                       /* FALLTHRU */
  case 2: ord_coef[2] = f0*(-1 + f2);                                 /* FALLTHRU */
  case 1: ord_coef[1] = 1 - f2;                                       /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                           break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(asin) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) < 1, "invalid domain asin("FMT")", VAL(a0));
  NUM f0 = asin(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // asin(x) = -i*ln(i*x + sqrt(1-x^2))
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_scl(c, -I, c);
#else
    ctpsa_t *t = GET_TMPR(c);
    mad_ctpsa_complex(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_scl(t, -I, t);
    mad_ctpsa_real(t, c);
    REL_TMPC(t);
#endif
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/sqrt(1-a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(5./16 + a2*(5./6 + 1./6*a2)) *f4*f4*f2*f1; /* FALLTHRU */
  case 5: ord_coef[5] = (3./40 + a2*(3./5 + 1./5*a2)) *f4*f4*f1;       /* FALLTHRU */
  case 4: ord_coef[4] = a0*(3./8 + 1./4*a2) *f4*f2*f1;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./6 + 1./3*a2) *f4*f1;                       /* FALLTHRU */
  case 2: ord_coef[2] = a0*(1./2) *f2*f1;                              /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                            /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                            break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(acos) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) < 1, "invalid domain acos("FMT")", VAL(a0));
  NUM f0 = acos(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acos(x) = -i*ln(x+i*sqrt(1-x^2)) = -asin(x)+pi/2
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_axpb(I, c, M_PI_2, c);
#else
    ctpsa_t *t = GET_TMPR(c);
    mad_ctpsa_complex(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_axpb(I, t, M_PI_2, t);
    mad_ctpsa_real(t, c);
    REL_TMPC(t);
#endif
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/sqrt(1-a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = -a0*(5./16 + a2*(5./6 + 1./6*a2)) *f4*f4*f2*f1; /* FALLTHRU */
  case 5: ord_coef[5] = -(3./40 + a2*(3./5 + 1./5*a2)) *f4*f4*f1;       /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(3./8 + 1./4*a2) *f4*f2*f1;                 /* FALLTHRU */
  case 3: ord_coef[3] = -(1./6 + 1./3*a2) *f4*f1;                       /* FALLTHRU */
  case 2: ord_coef[2] = -a0*(1./2) *f2*f1;                              /* FALLTHRU */
  case 1: ord_coef[1] = -f1;                                            /* FALLTHRU */
  case 0: ord_coef[0] =  f0;                                            break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(atan) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = atan(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // atan(x) = i/2 ln((i+x) / (i-x))
#ifdef MAD_CTPSA_IMPL
    ctpsa_t *tn = GET_TMPC(c), *td = GET_TMPC(c);
    mad_ctpsa_copy(a, tn);
    mad_ctpsa_axpb(-1, tn, I, td);
    mad_ctpsa_set0(tn, 1, I);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, I/2, c);
#else
    ctpsa_t *tn = GET_TMPR(c), *td = GET_TMPR(c);
    mad_ctpsa_complex(a, NULL, tn);
    mad_ctpsa_axpb(-1, tn, I, td);
    mad_ctpsa_set0(tn, 1, I);
    mad_ctpsa_logxdy(tn, td, tn);
    mad_ctpsa_scl(tn, I/2, tn);
    mad_ctpsa_real(tn, c);
#endif
    REL_TMPC(td), REL_TMPC(tn);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1+a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = -a0*(1 + a2*(-10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(-2 + a2)) *f4*f1;      /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(-1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (-1./3 + a2) *f2*f1;               /* FALLTHRU */
  case 2: ord_coef[2] = -a0 *f2;                           /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(acot) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain acot("FMT")", VAL(a0));
  NUM f0 = atan(1/a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acot(x) = i/2 ln((x-i) / (x+i))
#ifdef MAD_CTPSA_IMPL
    ctpsa_t *tn = GET_TMPC(c), *td = GET_TMPC(c);
    mad_ctpsa_copy( a, tn);
    mad_ctpsa_copy(tn, td);
    mad_ctpsa_set0(tn, 1, -I);
    mad_ctpsa_set0(td, 1,  I);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, I/2, c);
#else
    ctpsa_t *tn = GET_TMPR(c), *td = GET_TMPR(c);
    mad_ctpsa_complex(a, NULL, tn);
    mad_ctpsa_copy(tn, td);
    mad_ctpsa_set0(tn, 1, -I);
    mad_ctpsa_set0(td, 1,  I);
    mad_ctpsa_logxdy(tn, td, tn);
    mad_ctpsa_scl(tn, I/2, tn);
    mad_ctpsa_real(tn, c);
#endif
    REL_TMPC(td), REL_TMPC(tn);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1+a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(1 + a2*(-10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (-1./5 + a2*(2 - a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(-1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./3 - a2) *f2*f1;               /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                           /* FALLTHRU */
  case 1: ord_coef[1] = -f1;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;                               break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(asinh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = asinh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // asinh(x) = log(x + sqrt(x^2+1))
    FUN(logaxpsqrtbpcx2)(a, 1, 1, 1, c);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/sqrt(a2+1), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(-5./16 + a2*(5./6 - 1./6*a2)) *f4*f4*f2*f1; /* FALLTHRU */
  case 5: ord_coef[5] = (3./40 + a2*(-3./5 + 1./5*a2)) *f4*f4*f1;       /* FALLTHRU */
  case 4: ord_coef[4] = a0*(3./8 - 1./4*a2) *f4*f2*f1;                  /* FALLTHRU */
  case 3: ord_coef[3] = (-1./6 + 1./3*a2) *f4*f1;                       /* FALLTHRU */
  case 2: ord_coef[2] = a0*(-1./2) *f2*f1;                              /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                             /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                             break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(acosh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 1, 1), "invalid domain acosh("FMT")", VAL(a0));
  NUM f0 = acosh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acosh(x) = ln(x + sqrt(x^2-1))
    FUN(logaxpsqrtbpcx2)(a, 1, -1, 1, c);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/sqrt(a2-1), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = -a0*(5./16 + a2*(5./6 + 1./6*a2)) *f4*f4*f2*f1; /* FALLTHRU */
  case 5: ord_coef[5] = (3./40 + a2*(3./5 + 1./5*a2)) *f4*f4*f1;        /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(3./8 + 1./4*a2) *f4*f2*f1;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./6 + 1./3*a2) *f4*f1;                        /* FALLTHRU */
  case 2: ord_coef[2] = -a0*(1./2) *f2*f1;                              /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                             /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                             break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(atanh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(< 1, != 1), "invalid domain atanh("FMT")", VAL(a0));
  NUM f0 = atanh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // atanh(x) = 1/2 ln((1+x) / (1-x))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(set0)(tn, 1, 1);
    FUN(axpb)(-1, a, 1, td);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1-a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(1 + a2*(10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                          /* FALLTHRU */
  case 1: ord_coef[1] = f1;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;                              break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(acoth) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(> 1, != 1 && a0 != 0), "invalid domain acoth("FMT")", VAL(a0));
  NUM f0 = atanh(1/a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    // acoth(x) = 1/2 ln((x+1) / (x-1))
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

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1-a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(1 + a2*(10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                          /* FALLTHRU */
  case 1: ord_coef[1] = f1;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;                              break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

void
FUN(erf) (const T *a, T *c)
{
  // ERF(X) is the integral from 0 to x from [2/sqrt(pi) * exp(-x*x)]
  assert(a && c);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
#ifdef MAD_CTPSA_IMPL
  cnum_t f0 = mad_cnum_erf(a0, 0);
#else
  num_t  f0 = mad_num_erf (a0);
#endif

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) { FUN(scalar)(c,f0); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    ensure(0, "NYI");
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = M_2_SQRTPI/exp(a2);

  switch(to) {
  case 6: ord_coef[6] = -a0*(1./6 + a2*(-2./9 + 2./45*a2)) *f1; /* FALLTHRU */
  case 5: ord_coef[5] = (1./10 + a2*(-2./5 + 2./15*a2)) *f1;    /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(-1./2 + 1./3*a2) *f1;              /* FALLTHRU */
  case 3: ord_coef[3] = (-1./3 + 2./3*a2) *f1;                  /* FALLTHRU */
  case 2: ord_coef[2] = -a0 *f1;                                /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                     /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                     break;
  assert(!"unexpected missing coefficients");
  }

  fun_fix_point(a,c,to,ord_coef);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(inv_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(inv)(a, CNUM(v), c); }

void FUN(invsqrt_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(invsqrt)(a, CNUM(v), c); }

void FUN(pown_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(pown)(a, CNUM(v), c); }

#endif

// --- end --------------------------------------------------------------------o
