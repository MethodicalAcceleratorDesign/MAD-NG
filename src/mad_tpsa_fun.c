/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA functions module implementation
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

#ifdef MAD_TPSA_TAYLOR_HORNER

static inline void
fun_taylor (const T *a, T *c, ord_t n, const NUM ord_coef[n+1])
{
  assert(a && c && ord_coef);
  assert(n >= 1); // ord 0 treated outside

  T *acp = GET_TMPX(c);
  FUN(copy)(a,acp);                 // copy of a
  FUN(set0)(acp,0,0);               // (a-a_0)
  FUN(setvar)(c,ord_coef[n],0,0);   // f(a_n)

  // Honer's method (slower by 50% - 100% because mul is always full order)
  while (n-- > 0) {
    FUN(mul)(acp,c,c);            //                    f^(n)(a_n)*(a-a_0)
    FUN(set0)(c,1,ord_coef[n]);   // f^(n-1)(a_{n-1}) + f^(n)(a_n)*(a-a_0)
  }
  REL_TMPX(acp);
}

#else

static inline void
fun_taylor (const T *a, T *c, ord_t n, const NUM ord_coef[n+1])
{
  assert(a && c && ord_coef);
  assert(n >= 1); // ord 0 treated outside

  T *acp;
  if (n >= 2) acp = GET_TMPX(c), FUN(copy)(a,acp);

  // n=1
  FUN(scl)(a, ord_coef[1], c);
  FUN(set0)(c, 0, ord_coef[0]);    // f(a) + f'(a)(a-a0)

  // n=2
  if (n >= 2) {
    T *pow = GET_TMPX(c);
    FUN(set0)(acp,0,0);            //  a-a0
    FUN(mul)(acp,acp,pow);         // (a-a0)^2
    FUN(acc)(pow,ord_coef[2],c);   // f(a0) + f'(a0)(a-a0) + f"(a0)(a-a0)^2

    // i=3..n
    if (n >= 3) {
      T *tmp = GET_TMPX(c), *t;

      for (ord_t i = 3; i <= n; ++i) {
        FUN(mul)(acp,pow,tmp);
        FUN(acc)(tmp,ord_coef[i],c); // f(a0) + ... + f^(i)(a0)(a-a0)^i
        SWAP(pow,tmp,t);
      }

      if (n & 1) SWAP(pow,tmp,t); // enforce even number of swaps
      REL_TMPX(tmp);
    }
    REL_TMPX(pow), REL_TMPX(acp);
  }
}
#endif

static inline void
sincos_taylor (const T *a, T *s, T *c,
               ord_t n_s, const NUM sin_coef[n_s+1],
               ord_t n_c, const NUM cos_coef[n_c+1])
{
  assert(a && s && c && sin_coef && cos_coef);
  assert(n_s >= 1 && n_c >= 1);

  T *acp;
  ord_t n = MAX(n_s,n_c);
  if (n >= 2) acp = GET_TMPX(c), FUN(copy)(a,acp);

  // n=1
  FUN(scl)(a, sin_coef[1], s); FUN(set0)(s, 0, sin_coef[0]);
  FUN(scl)(a, cos_coef[1], c); FUN(set0)(c, 0, cos_coef[0]);

  // n=2
  if (n >= 2) {
    T *pow = GET_TMPX(c);
    FUN(set0)(acp,0,0);
    FUN(mul)(acp,acp,pow);
    if (n_s >= 2) FUN(acc)(pow,sin_coef[2],s);
    if (n_c >= 2) FUN(acc)(pow,cos_coef[2],c);

    // i=3..n
    if (n >= 3) {
      T *tmp = GET_TMPX(c), *t;

      for (ord_t i = 3; i <= n; ++i) {
        FUN(mul)(acp,pow,tmp);
        if (n_s >= i) FUN(acc)(tmp,sin_coef[i],s);
        if (n_c >= i) FUN(acc)(tmp,cos_coef[i],c);
        SWAP(pow,tmp,t);
      }

      if (n & 1) SWAP(pow,tmp,t); // enforce even number of swaps
      REL_TMPX(tmp);
    }
    REL_TMPX(pow), REL_TMPX(acp);
  }
}

// --- public -----------------------------------------------------------------o

void
FUN(taylor) (const T *a, ssz_t n, const NUM coef[n], T *c)
{
  assert(a && c && coef); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ensure(n > 0, "invalid number of coefficients (>0 expected)");

  ord_t to = MIN3(n-1,c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,coef[0],0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  fun_taylor(a,c,to,coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(inv) (const T *a, NUM v, T *c) // c = v/a    // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain inv("FMT")", VAL(a0));
  NUM f0 = 1/a0;

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,v*f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * f0;

  fun_taylor(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(invsqrt) (const T *a, NUM v, T *c) // v/sqrt(a),checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain invsqrt("FMT")", VAL(a0));
  NUM f0 = 1/sqrt(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,v*f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-1);

  fun_taylor(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sqrt) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain sqrt("FMT")", VAL(a0));
  NUM f0 = sqrt(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.0*o) * (2.0*o-3);

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(exp) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = exp(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (int o = 1; o <= to; ++o)
    ord_coef[o] = ord_coef[o-1] / o;

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(log) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain log("FMT")", VAL(a0));
  NUM f0 = log(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1], _a0 = 1/a0;
  ord_coef[0] = f0;
  ord_coef[1] = _a0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / o * (o-1);

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(pow) (const T *a, const T *b, T *c)          // checked for real and complex
{
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b);
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(mul)(b,t,c);
  FUN(exp)(c,c);
  REL_TMPX(t);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(pown) (const T *a, NUM v, T *c)              // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(scl)(t,v,c);
  FUN(exp)(c,c);
  REL_TMPX(t);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sincos) (const T *a, T *s, T *c)             // checked for real and complex
{
  assert(a && s && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == s->d && a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sin(a0), ca = cos(a0);

  if (a->hi == 0) {
    FUN(setvar)(s, sa, 0, 0);
    FUN(setvar)(c, ca, 0, 0);
    DBGTPSA(c); DBGTPSA(s); DBGFUN(<-);
    return;
  }

  ord_t sto = MIN(s->mo,s->d->to),
        cto = MIN(c->mo,c->d->to);
  if (!sto || !cto) {
    if (!sto) FUN(setvar)(s, sa, 0, 0);
    else      FUN(sin)(a,s);
    if (!cto) FUN(setvar)(c, ca, 0, 0);
    else      FUN(cos)(a,c);
    DBGTPSA(c); DBGTPSA(s); DBGFUN(<-);
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

  sincos_taylor(a,s,c, sto,sin_coef, cto,cos_coef);
  DBGTPSA(c); DBGTPSA(s); DBGFUN(<-);
}

void
FUN(sin) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sin(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cos(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(cos) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cos(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = -sin(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(tan) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(cos(a0) != 0, "invalid domain tan("FMT")", VAL(a0));
  NUM f0 = tan(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(cot) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(sin(a0) != 0, "invalid domain cot("FMT")", VAL(a0));
  NUM f0 = tan(M_PI_2 - a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  T *t = GET_TMPX(c);
  FUN(sincos)(a,t,c);
  FUN(div)(c,t,c);
  REL_TMPX(t);
  DBGTPSA(c); DBGFUN(<-);
  return;

#if 0
  // Inaccurate expansion for small a0
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c);
#endif
}

void
FUN(sinc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  ord_t to = MIN(c->mo,c->d->to);

  if (!to || a->hi == 0) {
#ifdef MAD_CTPSA_IMPL
    cnum_t f0 = mad_cnum_sinc(a0);
#else
    num_t  f0 = mad_num_sinc (a0);
#endif
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];

  if (fabs(a0) > 1e-12) { // sin(x)/x
    T *t = GET_TMPX(c);
    FUN(sin)(a,t);
    FUN(div)(t,a,c);
    REL_TMPX(t);
    DBGTPSA(c); DBGFUN(<-);
    return;
  }
// prefer explicit above? not better in term of stability...
//    NUM sa = sin(a0), ca = cos(a0), _a0 = 1/a0, f1;
//    num_t fo = 1;
//    ord_coef[0] = f0;
//    ord_coef[1] = (ca - f0)*_a0;
//    for (int o = 2; o <= to; ++o) {
//      fo *= o; // formula numerically unstable in (0, 0.5), need some work
//      f1  = o & 1 ? (ca=-ca,ca) : (sa=-sa,sa);
//      ord_coef[o] = (f1/fo - ord_coef[o-1])*_a0;
//      // printf("[%02d]=%+.17e\n", o, 1 - ord_coef[o-1]*fo/f1);
//    }

  // sinc(x)
  ord_coef[0] = 1;
  ord_coef[1] = 0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o * (o+1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sincosh) (const T *a, T *sh, T *ch)          // checked for real and complex
{
  assert(a && sh && ch); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == sh->d && a->d == ch->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);

  if (a->hi == 0) {
    FUN(setvar)(sh, sa, 0, 0);
    FUN(setvar)(ch, ca, 0, 0);
    DBGTPSA(ch); DBGTPSA(sh); DBGFUN(<-);
    return;
  }

  ord_t sto = MIN(sh->mo,sh->d->to),
        cto = MIN(ch->mo,ch->d->to);
  if (!sto || !cto) {
    if (!sto) FUN(setvar)(sh, sa, 0, 0);
    else      FUN(sinh)(a,sh);
    if (!cto) FUN(setvar)(ch, ca, 0, 0);
    else      FUN(cosh)(a,ch);
    DBGTPSA(ch); DBGTPSA(sh); DBGFUN(<-);
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

  sincos_taylor(a,sh,ch, sto,sin_coef, cto,cos_coef);
  DBGTPSA(ch); DBGTPSA(sh); DBGFUN(<-);
}

void
FUN(sinh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sinh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cosh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(cosh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cosh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = sinh(a0);
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(tanh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(coth) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);
  ensure(f0 != 0, "invalid domain coth("FMT")", VAL(a0));
  f0 = 1/f0;

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sinhc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  ord_t to = MIN(c->mo,c->d->to);

  if (!to || a->hi == 0) {
#ifdef MAD_CTPSA_IMPL
    cnum_t f0 = mad_cnum_sinhc(a0);
#else
    num_t  f0 = mad_num_sinhc (a0);
#endif
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];

  if (fabs(a0) > 1e-12) { // sinh(x)/x
    T *t = GET_TMPX(c);
    FUN(sinh)(a,t);
    FUN(div)(t,a,c);
    REL_TMPX(t);
    DBGTPSA(c); DBGFUN(<-);
    return;
  }

  // sinc(x)
  ord_coef[0] = 1;
  ord_coef[1] = 0;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o * (o+1));

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(asin) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) < 1, "invalid domain asin("FMT")", VAL(a0));
  NUM f0 = asin(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

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
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acos) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) < 1, "invalid domain acos("FMT")", VAL(a0));
  NUM f0 = acos(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

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
    DBGTPSA(c); DBGFUN(<-);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = -1/sqrt(1-a2), f2 = f1*f1, f4 = f2*f2;
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(atan) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = atan(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

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
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acot) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain acot("FMT")", VAL(a0));
  NUM f0 = atan(1/a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

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
    DBGTPSA(c); DBGFUN(<-);
    return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = -1/(1+a2), f2 = f1*f1, f4 = f2*f2;
  switch(to) {
  case 6: ord_coef[6] = a0*(1 + a2*(-10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(-2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(-1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (-1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                           /* FALLTHRU */
  case 1: ord_coef[1] = f1;                               /* FALLTHRU */
  case 0: ord_coef[0] = f0;                               break;
  assert(!"unexpected missing coefficients");
  }

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(asinh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = asinh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    // asinh(x) = log(x + sqrt(x^2+1))
    FUN(logaxpsqrtbpcx2)(a, 1, 1, 1, c);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acosh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 1, 1), "invalid domain acosh("FMT")", VAL(a0));
  NUM f0 = acosh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    // acosh(x) = ln(x + sqrt(x^2-1))
    FUN(logaxpsqrtbpcx2)(a, 1, -1, 1, c);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(atanh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(< 1, != 1), "invalid domain atanh("FMT")", VAL(a0));
  NUM f0 = atanh(a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  if (to > MANUAL_EXPANSION_ORD) {
    // atanh(x) = 1/2 ln((1+x) / (1-x))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(set0)(tn, 1, 1);
    FUN(axpb)(-1, a, 1, td);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn);
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acoth) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(> 1, != 1 && a0 != 0), "invalid domain acoth("FMT")", VAL(a0));
  NUM f0 = atanh(1/a0);

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

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
    DBGTPSA(c); DBGFUN(<-);
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

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(erf) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  ensure(a->d == c->d, "incompatible GTPSA (descriptors differ)");

  // erf(z) = 2/sqrt(pi) \int_0^z exp(-t^2) dt
  NUM a0 = a->coef[0];
#ifdef MAD_CTPSA_IMPL
  cnum_t f0 = mad_cnum_erf(a0, 0);
#else
  num_t  f0 = mad_num_erf (a0);
#endif

  ord_t to = MIN(c->mo,c->d->to);
  if (!to || a->hi == 0) {
    FUN(setvar)(c,f0,0,0); DBGTPSA(c); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = M_2_SQRTPI*exp(-a2);
  ord_coef[0] = f0;
  ord_coef[1] = f1;
  for (int o = 2; o <= to; ++o)
    ord_coef[o] = -2*((o-2)*ord_coef[o-2]/(o-1) + ord_coef[o-1]*a0) / o;

  fun_taylor(a,c,to,ord_coef);
  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(erfc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a);
  FUN(erf)(a,c);
  FUN(axpb)(-1,c,1,c);
  DBGTPSA(c); DBGFUN(<-);
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

/*
Recurrence of Taylor coefficients:
----------------------------------
(f(g(x)))'   = f'(g(x)).g'(x)
(f(x).g(x))' = f'(x)g(x) + f(x)g'(x)

-- sinc(z) -----
[0]  sinc(z)      =   sin(z)/z ; cos(z)/z = [c] ; sin(z) = sz ; cos(z) = cz
[1] (sinc(z))'    =  cz/z -1sz/z^2                                                     =  [c]-1*[0]/z
[2] (sinc(z))''   = -sz/z -2cz/z^2 +2!sz/z^3                                           = -[0]-2*[1]/z
[3] (sinc(z))'''  = -cz/z +3sz/z^2 +3!cz/z^3 - 3!sz/z^4                                = -[c]-3*[2]/z
[4] (sinc(z))^(4) =  sz/z +4cz/z^2 -12sz/z^3 - 4!cz/z^4 + 4!sz/z^5                     =  [0]-4*[3]/z
[5] (sinc(z))^(5) =  cz/z -5sz/z^2 -20cz/z^3 + 60sz/z^4 + 5!cz/z^5 -5!sz/z^6           =  [c]-5*[4]/z
[6] (sinc(z))^(6) = -sz/z -6cz/z^2 +30sz/z^3 +120cz/z^4 -360sz/z^5 -6!cz/z^6 +6!sz/z^7 = -[0]-6*[5]/z

-- erf(z) -----
[0]  erf(z)
[1] (erf(z))'    =     1
[2] (erf(z))''   = -2*      z                                 = -2*(0*[0]+[1]*z)
[3] (erf(z))'''  = -2*(1 -2*z^2)                              = -2*(1*[1]+[2]*z)
[4] (erf(z))^(4) = -2*(-4*z -2*(1-2*z^2)*z)                   = -2*(2*[2]+[3]*z)
[5] (erf(z))^(5) = -2*(-6*(1-2*z^2) +4*(3*z-2*z^3)*z)         = -2*(3*[3]+[4]*z)
[6] (erf(z))^(6) = -2*(16*(3*z-2*z^3) +4*(3-12*z^2+4*z^4)*z)  = -2*(4*[4]+[5]*z)
                   % *exp(-z^2) *2/sqrt(pi)
{0} = 0
{1} = 1 *exp(-z^2)
(exp(-z^2))'
    = exp'(-z^2).(-z^2)'
{2} = -2*z *exp(-z^2)                                         = -2*(0*{0}+{1}*z)
-2*(z*exp(-z^2))'
    = -2*(z'*exp(-z^2) + z*(exp(-z^2))')
    = -2*(exp(-z^2) + z*(-2*z*exp(-z^2)))
{3} = -2* (1 -2*z^2) *exp(-z^2)                               = -2*(1*{1}+{2}*z)
-2*((1-2*z^2)*exp(-z^2))' =
    = -2*((1-2*z^2)'*exp(-z^2) + (1-2*z^2)*(exp(-z^2))')
    = -2*(-4*z*exp(-z^2) + (1-2*z^2)*(-2*z*exp(-z^2)))
{4} = -2*(-4*z -2*(1-2*z^2)*z) *exp(-z^2)                     = -2*(2*{2}+{3}*z)
    =  4*(3*z-2*z^3) *exp(-z^2)
4*((3*z-2*z^3)*exp(-z^2))' =
    = 4*((3*z-2*z^3)'*exp(-z^2) + (3*z-2*z^3)*(exp(-z^2))')
    = 4*(3*(1-2*z^2)*exp(-z^2) + (3*z-2*z^3)*(-2*z*exp(-z^2)))
{5} = -2*(3*-2*(1-2*z^2) + 4*(3*z-2*z^3)*z) *exp(-z^2)        = -2*(3*{3}+{4}*z)
    = 4*(3-12*z^2+4*z^4) *exp(-z^2)
4*((3-12*z^2+4*z^4)*exp(-z^2))' =
    = 4*((3-12*z^2+4*z^4)'*exp(-z^2) + (3-12*z^2+4*z^4)*(exp(-z^2))')
    = 4*((-24*z+16*z^3)*exp(-z^2) + (3-12*z^2+4*z^4)*(-2*z*exp(-z^2)))
    = 4*(-2*(12*z-8*z^3) -2*(3-12*z^2+4*z^4)*z) *exp(-z^2)
{6} = -2*(16*(3*z-2*z^3) +4*(3-12*z^2+4*z^4)*z) *exp(-z^2)    = -2*(4*{4}+{5}*z)
    = 4*(-30*z+40*z^3-8*z^5) *exp(-z^2)
...
*/
