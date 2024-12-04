/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA functions module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: A. E. Scurria (Numerical stability for high orders)
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o
*/

#include <string.h>
#include <float.h>
#include <complex.h>

#include "mad_log.h"
#include "mad_cst.h"
#include "mad_num.h"
#include "mad_tpsa_impl.h"
#include "mad_ctpsa_impl.h"

// --- local ------------------------------------------------------------------o

#define OLD_SERIES 0
#define DBG_SERIES 0

enum { MANUAL_EXPANSION_ORD = 6 };

static inline void
fun_taylor_h (const T *a, T *c, ord_t n, const NUM ord_coef[n+1])
{
  assert(a && c && ord_coef);
  assert(n >= 1); // ord 0 treated outside

  T *acp = GET_TMPX(c);
  FUN(copy)(a,acp);               // copy of a
  FUN(seti)(acp,0,0,0);           // (a-a_0)
  FUN(setval)(c,ord_coef[n]);     // f(a_n)

  // Horner's method (slower by 50% - 100% because mul is always full order)
  while (n-- > 0) {
    FUN(mul)(acp,c,c);            //                    f^(n)(a_n)*(a-a_0)
    FUN(seti)(c,0,1,ord_coef[n]); // f^(n-1)(a_{n-1}) + f^(n)(a_n)*(a-a_0)
  }
  REL_TMPX(acp);
}

static inline void
fun_taylor (const T *a, T *c, ord_t n, const NUM ord_coef[n+1])
{
  assert(a && c && ord_coef);
  assert(n >= 1); // ord 0 treated outside

  T *acp;
  if (n >= 2) acp = GET_TMPX(c), FUN(copy)(a,acp);

  // n=1
  FUN(scl)(a, ord_coef[1], c);
  FUN(seti)(c, 0, 0, ord_coef[0]); // f(a0) + f'(a0)(a-a0)

  // n=2
  if (n >= 2) {
    T *pow = GET_TMPX(c);
    FUN(seti)(acp,0,0,0);          //  a-a0
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

static inline void
sincos_taylor (const T *a, T *s, T *c,
               ord_t n_s, const NUM sin_coef[n_s+1],
               ord_t n_c, const NUM cos_coef[n_c+1])
{
  assert(a && s && c && sin_coef && cos_coef);
  assert(n_s >= 1 && n_c >= 1);

  ord_t n = MAX(n_s,n_c);
  T *acp = GET_TMPX(c); FUN(copy)(a,acp);

  // n=1
  FUN(scl)(acp, sin_coef[1], s); FUN(seti)(s, 0, 0, sin_coef[0]);
  FUN(scl)(acp, cos_coef[1], c); FUN(seti)(c, 0, 0, cos_coef[0]);

  // n=2
  if (n >= 2) {
    T *pow = GET_TMPX(c);
    FUN(seti)(acp,0,0,0);
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
    REL_TMPX(pow);
  }
  REL_TMPX(acp);
}

// --- public -----------------------------------------------------------------o

void
FUN(taylor_h) (const T *a, ssz_t n, const NUM coef[n], T *c)
{
  assert(a && c && coef); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  ensure(n > 0, "invalid number of coefficients (>0 expected)");

  ord_t to = MIN(n-1, c->mo);
  if (!to || FUN(isval)(a)) { FUN(setval)(c,coef[0]); DBGFUN(<-); return; }

  fun_taylor_h(a,c,to,coef);
  DBGFUN(<-);
}

void
FUN(taylor) (const T *a, ssz_t n, const NUM coef[n], T *c)
{
  assert(a && c && coef); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  ensure(n > 0, "invalid number of coefficients (>0 expected)");

  ord_t to = MIN(n-1, c->mo);
  if (!to || FUN(isval)(a)) { FUN(setval)(c,coef[0]); DBGFUN(<-); return; }

  fun_taylor(a,c,to,coef);
  DBGFUN(<-);
}

void
FUN(inv) (const T *a, NUM v, T *c) // c = v/a    // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain inv("FMT")", VAL(a0));

  NUM f0 = NUMF(inv)(a0);
  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,v*f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * f0;

  fun_taylor(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
  DBGFUN(<-);
}

void
FUN(invsqrt) (const T *a, NUM v, T *c) // v/sqrt(a),checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain invsqrt("FMT")", VAL(a0));

  NUM _a0 = NUMF(inv)(a0);
  NUM  f0 = NUMF(inv)(sqrt(a0));
  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,v*f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.*o) * (2.*o-1);

  fun_taylor(a,c,to,ord_coef);
  if (v != 1) FUN(scl)(c,v,c);
  DBGFUN(<-);
}

void
FUN(sqrt) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain sqrt("FMT")", VAL(a0));
  NUM f0 = sqrt(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM _a0 = NUMF(inv)(a0);
  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / (2.*o) * (2.*o-3);

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(exp) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM f0 = exp(a->coef[0]);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = ord_coef[o-1] / o;

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(log) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 0, a0 != 0), "invalid domain log("FMT")", VAL(a0));
  NUM f0 = log(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM _a0 = NUMF(inv)(a0);
  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = _a0;
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-1] * _a0 / o * (o-1.);

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(pow) (const T *a, const T *b, T *c)          // checked for real and complex
{
  assert(a && b && c); DBGFUN(->);
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(mul)(b,t,c);
  FUN(exp)(c,c);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(pown) (const T *a, NUM v, T *c)              // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  T *t = GET_TMPX(c);
  FUN(log)(a,t);
  FUN(scl)(t,v,c);
  FUN(exp)(c,c);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(sincos) (const T *a, T *s, T *c)             // checked for real and complex
{
  assert(a && s && c); DBGFUN(->);
  ensure(IS_COMPAT(a,s,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sin(a0), ca = cos(a0);

  if (a->hi == 0) {
    FUN(setval)(s, sa);
    FUN(setval)(c, ca);
    DBGFUN(<-); return;
  }

  ord_t sto = s->mo, cto = c->mo;
  if (!sto || !cto) {
    if (!sto) FUN(setval)(s, sa);
    else      FUN(sin)(a,s);
    if (!cto) FUN(setval)(c, ca);
    else      FUN(cos)(a,c);
    DBGFUN(<-); return;
  }

  // ord 0, 1
  NUM sin_coef[sto+1], cos_coef[cto+1];
  sin_coef[0] = sa;  cos_coef[0] =  ca;
  sin_coef[1] = ca;  cos_coef[1] = -sa;

  // ords 2..to
  for (ord_t o = 2; o <= sto; ++o )
    sin_coef[o] = -sin_coef[o-2] / (o*(o-1.));
  for (ord_t o = 2; o <= cto; ++o )
    cos_coef[o] = -cos_coef[o-2] / (o*(o-1.));

  sincos_taylor(a,s,c, sto,sin_coef, cto,cos_coef);
  DBGFUN(<-);
}

void
FUN(sin) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sin(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cos(a0);
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1.));

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(cos) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cos(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = -sin(a0);
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = -ord_coef[o-2] / (o*(o-1.));

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(tan) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(cos(a0) != 0, "invalid domain tan("FMT")", VAL(a0));
  NUM f0 = tan(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t); DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(cot) (const T *a, T *c)                      // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(sin(a0) != 0, "invalid domain cot("FMT")", VAL(a0));
  NUM f0 = tan(M_PI_2 - a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincos)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t); DBGFUN(<-); return;
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

  fun_taylor(a,c,to,ord_coef);
}

void
FUN(sinc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  NUM f0 = NUMF(sinc)(a0);
  ord_t to = c->mo;

  if (!to || FUN(isval)(a)) {
    FUN(setval)(c,f0); DBGFUN(<-); return;
  }

  if (fabs(a0) > 0.75) { // sinc(x), |x| > 0.75
    T *t = GET_TMPX(c);
    FUN(sin)(a,t);
    FUN(div)(t,a,c);
    FUN(seti)(c,0,0,f0); // scalar part not very accurate
    REL_TMPX(t); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];

  if (fabs(a0) > 1e-12) { // sinc(x), 1e-12 < |x| <= 0.75
    num_t odd_coef[7] = {
      -3, 30, -840, 45360, -3991680, 518918400, -93405312000
    };
    int s = 1; // o/2: even = 1, odd = -1
    num_t fact = 1;
    for (ord_t o = 1; o <= to; o += 2, s = -s) {
      NUM v1 = 0, v2 = 0;
      FOR(i,7) {
        v1 += 1 / odd_coef[i] * NUMF(powi)(a0,2*i+1);
        v2 += 1 / odd_coef[i] * NUMF(powi)(a0,2*i)*(2*i+1);
      }
      fact *= o*(o+1.);
      ord_coef[ o          ] = s * v1 / fact*(o+1);
      ord_coef[(o+1)%(to+1)] = s * v2 / fact;
      if (o+2 <= to) {
        static const num_t stp_coef[7] = {
          -2, 12, -240, 10080, -725760, 79833600, -12454041600
        };
        FOR(i,7) odd_coef[i] += stp_coef[i];
      }
    }
    ord_coef[0] = f0;
  } else { // sinc(x), |x| <= 1e-12
    ord_coef[0] = 1;
    ord_coef[1] = 0;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = -ord_coef[o-2] / (o*(o+1.));
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(sincosh) (const T *a, T *sh, T *ch)          // checked for real and complex
{
  assert(a && sh && ch); DBGFUN(->);
  ensure(IS_COMPAT(a,sh,ch), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], sa = sinh(a0), ca = cosh(a0);

  if (a->hi == 0) {
    FUN(setval)(sh, sa);
    FUN(setval)(ch, ca);
    DBGFUN(<-); return;
  }

  ord_t sto = sh->mo, cto = ch->mo;
  if (!sto || !cto) {
    if (!sto) FUN(setval)(sh, sa);
    else      FUN(sinh)(a,sh);
    if (!cto) FUN(setval)(ch, ca);
    else      FUN(cosh)(a,ch);
    DBGFUN(<-); return;
  }

  // ord 0, 1
  NUM sin_coef[sto+1], cos_coef[cto+1];
  sin_coef[0] = sa;  cos_coef[0] = ca;
  sin_coef[1] = ca;  cos_coef[1] = sa;

  // ords 2..to
  for (ord_t o = 2; o <= sto; ++o)
    sin_coef[o] = sin_coef[o-2] / (o*(o-1.));
  for (ord_t o = 2; o <= cto; ++o)
    cos_coef[o] = cos_coef[o-2] / (o*(o-1.));

  sincos_taylor(a,sh,ch, sto,sin_coef, cto,cos_coef);
  DBGFUN(<-);
}

void
FUN(sinh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = sinh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = cosh(a0);
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1.));

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(cosh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = cosh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = sinh(a0);
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = ord_coef[o-2] / (o*(o-1.));

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(tanh) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(t,c,c);
    REL_TMPX(t); DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(coth) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = tanh(a0);
  ensure(f0 != 0, "invalid domain coth("FMT")", VAL(a0));

  f0 = NUMF(inv)(f0);
  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) {
    T *t = GET_TMPX(c);
    FUN(sincosh)(a,t,c);
    FUN(div)(c,t,c);
    REL_TMPX(t); DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(sinhc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  NUM f0 = NUMF(sinhc)(a0);
  ord_t to = c->mo;

  if (!to || FUN(isval)(a)) {
    FUN(setval)(c,f0); DBGFUN(<-); return;
  }

  if (fabs(a0) > 0.75) { // sinhc(x), |x| > 0.75
    T *t = GET_TMPX(c);
    FUN(sinh)(a,t);
    FUN(div)(t,a,c);
    FUN(seti)(c,0,0,f0); // scalar part not very accurate
    REL_TMPX(t); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];
  if (fabs(a0) > 1e-12) { // sinhc(x), 1e-12 < |x| <= 0.75
    num_t odd_coef[7] = {
      3, 30, 840, 45360, 3991680, 518918400, 93405312000
    };
    num_t fact = 1;
    for (ord_t o = 1; o <= to; o += 2) {
      NUM v1 = 0, v2 = 0;
      FOR(i,7) {
        v1 += 1 / odd_coef[i] * NUMF(powi)(a0,2*i+1);
        v2 += 1 / odd_coef[i] * NUMF(powi)(a0,2*i)*(2*i+1);
      }
      fact *= o*(o+1.);
      ord_coef[ o          ] = v1 / fact*(o+1);
      ord_coef[(o+1)%(to+1)] = v2 / fact;
      if (o+2 <= to) {
        static const num_t stp_coef[7] = {
          2, 12, 240, 10080, 725760, 79833600, 12454041600
        };
        FOR(i,7) odd_coef[i] += stp_coef[i];
      }
    }
    ord_coef[0] = f0;
  } else { // sinhc(x), |x| <= 1e-12
    ord_coef[0] = 1;
    ord_coef[1] = 0;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = ord_coef[o-2] / (o*(o+1.));
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(asin) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(fabs(a0) < 1, 1), "invalid domain asin("FMT")", VAL(a0));
  NUM f0 = asin(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // asin(x) = -i*ln(i*x + sqrt(1-x^2))
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_scl(c, -I, c);
#else
    ctpsa_t *t = GET_TMPC(c);
    mad_ctpsa_cplx(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_scl(t, -I, t);
    mad_ctpsa_real(t, c);
    REL_TMPC(t);
#endif
    DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(acos) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(fabs(a0) < 1, 1), "invalid domain acos("FMT")", VAL(a0));
  NUM f0 = acos(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) {  // use simpler and faster approach?
    // acos(x) = -i*ln(x + i*sqrt(1-x^2)) = -asin(x) + pi/2
#ifdef MAD_CTPSA_IMPL
    mad_ctpsa_logaxpsqrtbpcx2(a, I, 1, -1, c);
    mad_ctpsa_axpb(I, c, M_PI_2, c);
#else
    ctpsa_t *t = GET_TMPC(c);
    mad_ctpsa_cplx(a, NULL, t);
    mad_ctpsa_logaxpsqrtbpcx2(t, I, 1, -1, t);
    mad_ctpsa_axpb(I, t, M_PI_2, t);
    mad_ctpsa_real(t, c);
    REL_TMPC(t);
#endif
    DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(atan) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = atan(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

#if OLD_SERIES
  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // atan(x) = -i/2 ln((i-x) / (i+x)) ; pole = +/- 1i
#ifdef MAD_CTPSA_IMPL
    ctpsa_t *tn = GET_TMPX(c), *td = GET_TMPX(c);
    mad_ctpsa_axpb(-1, a, I, tn);
    mad_ctpsa_axpb( 1, a, I, td);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, -I/2, c);
#else
    ctpsa_t *tn = GET_TMPC(c), *td = GET_TMPC(c);
    mad_ctpsa_cplx(a, NULL, td);
    mad_ctpsa_axpb(-1, td, I, tn);
    mad_ctpsa_seti(td, 0, 1, I);
    mad_ctpsa_logxdy(tn, td, tn);
    mad_ctpsa_scl(tn, -I/2, tn);
    mad_ctpsa_real(tn, c);
#endif
    REL_TMPC(td), REL_TMPC(tn); DBGFUN(<-); return;
  }
#endif // OLD_SERIES

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1+a2), f2 = f1*f1, f4 = f2*f2;
  // fill orders <= MANUAL_EXPANSION_ORD
  switch(MIN(to, MANUAL_EXPANSION_ORD)) {
  case 6: ord_coef[6] = -a0*(1 + a2*(-10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(-2 + a2)) *f4*f1;      /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(-1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (-1./3 + a2) *f2*f1;               /* FALLTHRU */
  case 2: ord_coef[2] = -a0 *f2;                           /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                break;
  assert(!"unexpected missing coefficients");
  }
  // fill orders > MANUAL_EXPANSION_ORD
  for (ord_t o = 1+MANUAL_EXPANSION_ORD; o <= to; ++o) {
    ord_t n = (o-3)/2+1;
    int s = 2*(o & 1)-1; // o: even = -1, odd = 1
    NUM v = 0;
    for (ord_t i = 0; i <= n; ++i, s = -s) {
      num_t c = s * mad_num_powi(2,o-2*i-1) * mad_num_binom(o-i-1,i);
      v += c * NUMF(div)(NUMF(powi)(a0,o-2*i-1), NUMF(powi)(1+a2,o-i));
    }
    ord_coef[o] = v/o;
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(acot) (const T *a, T *c)                     // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(a0 != 0, "invalid domain acot("FMT")", VAL(a0));

  NUM f0 = atan(NUMF(inv)(a0));
  ord_t to = c->mo;

#if OLD_SERIES
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // acot(x) = i/2 ln((x-i) / (x+i))
#ifdef MAD_CTPSA_IMPL
    ctpsa_t *tn = GET_TMPX(c), *td = GET_TMPX(c);
    mad_ctpsa_copy( a, tn);
    mad_ctpsa_copy(tn, td);
    mad_ctpsa_seti(tn, 0, 1, -I);
    mad_ctpsa_seti(td, 0, 1,  I);
    mad_ctpsa_logxdy(tn, td, c);
    mad_ctpsa_scl(c, I/2, c);
#else
    ctpsa_t *tn = GET_TMPC(c), *td = GET_TMPC(c);
    mad_ctpsa_cplx(a, NULL, tn);
    mad_ctpsa_copy(tn, td);
    mad_ctpsa_seti(tn, 0, 1, -I);
    mad_ctpsa_seti(td, 0, 1,  I);
    mad_ctpsa_logxdy(tn, td, tn);
    mad_ctpsa_scl(tn, I/2, tn);
    mad_ctpsa_real(tn, c);
#endif
    REL_TMPC(td), REL_TMPC(tn); DBGFUN(<-); return;
  }
#endif // OLD_SERIES

  NUM ord_coef[to+1], a2 = a0*a0, f1 = -1/(1+a2), f2 = f1*f1, f4 = f2*f2;
  // fill orders <= MANUAL_EXPANSION_ORD
  switch(MIN(to, MANUAL_EXPANSION_ORD)) {
  case 6: ord_coef[6] = a0*(1 + a2*(-10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(-2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(-1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (-1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                           /* FALLTHRU */
  case 1: ord_coef[1] = f1;                               /* FALLTHRU */
  case 0: ord_coef[0] = f0;                               break;
  assert(!"unexpected missing coefficients");
  }
  // fill orders > MANUAL_EXPANSION_ORD
  for (ord_t o = 1+MANUAL_EXPANSION_ORD; o <= to; ++o) {
    ord_t n = (o-3)/2+1;
    int s = 2*(o & 1)-1; // o: even = -1, odd = 1
    NUM v = 0;
    for (ord_t i = 0; i <= n; ++i, s = -s) {
      num_t c = s * mad_num_powi(2,o-2*i-1) * mad_num_binom(o-i-1,i);
      v += c * NUMF(div)(NUMF(powi)(a0,o-2*i-1), NUMF(powi)(1+a2,o-i));
    }
    ord_coef[o] = -v/o;
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(asinc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  NUM f0 = NUMF(asinc)(a0);
  ord_t to = c->mo;

  if (!to || FUN(isval)(a)) {
    FUN(setval)(c,f0); DBGFUN(<-); return;
  }

  if (fabs(a0) > 0.58) { // asinc(x), |x| > 0.58
    T *t = GET_TMPX(c);
    FUN(asin)(a,t);
    FUN(div)(t,a,c);
    FUN(seti)(c,0,0,f0); // scalar part not very accurate
    REL_TMPX(t); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];

  if (fabs(a0) > 1e-12) { // asinc(x), 1e-12 < |x| <= 0.58
    FOR(i,to+1) ord_coef[i] = 0;

    enum { ord = DESC_MAX_ORD+1 };
    static num_t scl_coef[ord+1] = {0};
    if (!scl_coef[0]) {
      scl_coef[0] = 1./3;
      FOR(i,1,ord+1) scl_coef[i] = scl_coef[i-1] * SQR(2*i+1.)/(i*(4*i+6.));
    }

    num_t fact = 1;
    num_t tmp_coef[ord+1]; FOR(i,ord+1) tmp_coef[i] = scl_coef[i];
    for (ord_t o = 1; o <= to; o += 2) {
      NUM a0pi = 1, dc1, dc2;
      fact *= o*(o+1.);
      FOR(i,ord+1) {
        tmp_coef[i] *= o > 1 ? CUB(2.*i+o) / (2*i+o+2) : 1;
        ord_coef [ o          ] += (dc1=a0*SQR(a0pi)*tmp_coef[i]*(  o+1) / fact);
        ord_coef [(o+1)%(to+1)] += (dc2=   SQR(a0pi)*tmp_coef[i]*(2*i+1) / fact);
        if (fabs(dc1)/fabs(ord_coef[o])            < SQR(DBL_EPSILON) &&
            fabs(dc2)/fabs(ord_coef[(o+1)%(to+1)]) < SQR(DBL_EPSILON)) {
          if (DBG_SERIES)
            printf("asinc: i=%d, o=%d, to=%d, |a0|=%.16e, |dc|=%.16e\n", i, o, to, fabs(a0), MAX(fabs(dc1),fabs(dc2)));
          break;
        } else if (i == ord) {
          if (DBG_SERIES)
            printf("asinc: i=%d, o=%d, to=%d, |a0|=%.16e, |dc|=%.16e\n", i, o, to, fabs(a0), MAX(fabs(dc1),fabs(dc2)));
          trace(1, "asinc did not converge after %d iterations", i);
        }
        a0pi *= a0;
      }
    }
    ord_coef[0] = f0;
  } else { // asinc(x), |x| <= 1e-12
    ord_coef[0] = 1;
    ord_coef[1] = 0;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = (ord_coef[o-2] * SQR(o-1.)) / (o*(o+1.));
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(asinh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0], f0 = asinh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // asinh(x) = ln(x + sqrt(x^2+1))
    FUN(logaxpsqrtbpcx2)(a, 1, 1, 1, c);
    DBGFUN(<-); return;
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
  DBGFUN(<-);
}

void
FUN(acosh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(SELECT(a0 > 1, 1), "invalid domain acosh("FMT")", VAL(a0));
  NUM f0 = acosh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

#if OLD_SERIES
  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // acosh(x) = ln(x + sqrt(x^2-1))
    FUN(logaxpsqrtbpcx2)(a, 1, -1, 1, c);
    DBGFUN(<-); return;
  }
#endif // OLD_SERIES

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/sqrt(a2-1), f2 = f1*f1, f4 = f2*f2;
  // fill orders <= MANUAL_EXPANSION_ORD
  switch(MIN(to, MANUAL_EXPANSION_ORD)) {
  case 6: ord_coef[6] = -a0*(5./16 + a2*(5./6 + 1./6*a2)) *f4*f4*f2*f1; /* FALLTHRU */
  case 5: ord_coef[5] = (3./40 + a2*(3./5 + 1./5*a2)) *f4*f4*f1;        /* FALLTHRU */
  case 4: ord_coef[4] = -a0*(3./8 + 1./4*a2) *f4*f2*f1;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./6 + 1./3*a2) *f4*f1;                        /* FALLTHRU */
  case 2: ord_coef[2] = -a0*(1./2) *f2*f1;                              /* FALLTHRU */
  case 1: ord_coef[1] = f1;                                             /* FALLTHRU */
  case 0: ord_coef[0] = f0;                                             break;
  assert(!"unexpected missing coefficients");
  }
  // fill orders > MANUAL_EXPANSION_ORD
  if (to > MANUAL_EXPANSION_ORD) {
    NUM _a0p1 = 1/(a0+1);
    NUM _a0m1 = 1/(a0-1);
    NUM   den = NUMF(powi)(-2,MANUAL_EXPANSION_ORD-1) * sqrt((a0+1)*(a0-1));

    for (ord_t o = 1+MANUAL_EXPANSION_ORD; o <= to; ++o) {
      int n = (o-1)/2+1;
      NUM num = 0; den *= -2;
      FOR(i,n) {
        int d = o-2*i-1 ? 1 : 2;
        num_t c = (i ? mad_num_fact2(2*i-1) : 1.) / d;
        num += c * mad_num_fact2(2*o-2*i-3) * mad_num_binom(o-1,i) *
              (NUMF(powi)(_a0p1,o-i-1) * NUMF(powi)(_a0m1,i) +
               NUMF(powi)(_a0m1,o-i-1) * NUMF(powi)(_a0p1,i) );
      }
      ord_coef[o] = num/den/mad_num_fact(o);
    }
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(atanh) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(< 1, != 1), "invalid domain atanh("FMT")", VAL(a0));
  NUM f0 = atanh(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

#if OLD_SERIES
  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // atanh(x) = 1/2 ln((1+x) / (1-x))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(seti)(tn, 0, 1, 1);
    FUN(axpb)(-1, a, 1, td);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn); DBGFUN(<-); return;
  }
#endif // OLD_SERIES

  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1-a2), f2 = f1*f1, f4 = f2*f2;
  // fill orders <= MANUAL_EXPANSION_ORD
  switch(MIN(to, MANUAL_EXPANSION_ORD)) {
  case 6: ord_coef[6] = a0*(1 + a2*(10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                          /* FALLTHRU */
  case 1: ord_coef[1] = f1;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;                              break;
  assert(!"unexpected missing coefficients");
  }
  // fill orders > MANUAL_EXPANSION_ORD
  for (ord_t o = 1+MANUAL_EXPANSION_ORD; o <= to; ++o) {
    ord_t n = (o-3)/2+1;
    int s = 2*(o & 1)-1; // o: even = -1, odd = 1
    NUM v = 0;
    for (ord_t i = 0; i <= n; ++i, s = -s) {
      num_t c = s * mad_num_powi(2,o-2*i-1) * mad_num_binom(o-i-1,i);
      v += c * NUMF(div)(NUMF(powi)(a0,o-2*i-1), NUMF(powi)(a2-1,o-i));
    }
    ord_coef[o] = -v/o;
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(acoth) (const T *a, T *c)                    // checked for real and complex
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  ensure(fabs(a0) SELECT(> 1, != 1 && a0 != 0), "invalid domain acoth("FMT")", VAL(a0));

  NUM f0 = atanh(NUMF(inv)(a0));
  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

#if OLD_SERIES
  if (to > MANUAL_EXPANSION_ORD) { // use simpler and faster approach?
    // acoth(x) = 1/2 ln((x+1) / (x-1))
    T *tn = GET_TMPX(c), *td = GET_TMPX(c);
    FUN(copy)(a, tn);
    FUN(seti)(tn, 0, 1, 1);
    FUN(copy)(a, td);
    FUN(seti)(td, 0, 1, -1);
    FUN(logxdy)(tn, td, c);
    FUN(scl)(c, 0.5, c);
    REL_TMPX(td), REL_TMPX(tn); DBGFUN(<-); return;
  }
#endif // OLD_SERIES
  
  NUM ord_coef[to+1], a2 = a0*a0, f1 = 1/(1-a2), f2 = f1*f1, f4 = f2*f2;
  // fill orders <= MANUAL_EXPANSION_ORD
  switch(MIN(to, MANUAL_EXPANSION_ORD)) {
  case 6: ord_coef[6] = a0*(1 + a2*(10./3 + a2)) *f4*f2; /* FALLTHRU */
  case 5: ord_coef[5] = (1./5 + a2*(2 + a2)) *f4*f1;     /* FALLTHRU */
  case 4: ord_coef[4] = a0*(1 + a2) *f4;                 /* FALLTHRU */
  case 3: ord_coef[3] = (1./3 + a2) *f2*f1;              /* FALLTHRU */
  case 2: ord_coef[2] = a0 *f2;                          /* FALLTHRU */
  case 1: ord_coef[1] = f1;                              /* FALLTHRU */
  case 0: ord_coef[0] = f0;                              break;
  assert(!"unexpected missing coefficients");
  }
  // fill orders > MANUAL_EXPANSION_ORD
  for (ord_t o = 1+MANUAL_EXPANSION_ORD; o <= to; ++o) {
    ord_t n = (o-3)/2+1;
    int s = 2*(o & 1)-1; // o: even = -1, odd = 1
    NUM v = 0;
    for (ord_t i = 0; i <= n; ++i, s = -s) {
      num_t c = s * mad_num_powi(2,o-2*i-1) * mad_num_binom(o-i-1,i);
      v += c * NUMF(div)(NUMF(powi)(a0,o-2*i-1), NUMF(powi)(a2-1,o-i));
    }
    ord_coef[o] = -v/o;
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(asinhc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  NUM f0 = NUMF(asinhc)(a0);
  ord_t to = c->mo;

  if (!to || FUN(isval)(a)) {
    FUN(setval)(c,f0); DBGFUN(<-); return;
  }

  if (fabs(a0) > 0.62) { // asinhc(x), |x| > 0.62
    T *t = GET_TMPX(c);
    FUN(asinh)(a,t);
    FUN(div)(t,a,c);
    FUN(seti)(c,0,0,f0); // scalar part not very accurate
    REL_TMPX(t); DBGFUN(<-); return;
  }

  NUM ord_coef[to+1];

  if (fabs(a0) > 1e-12) { // asinhc(x), 1e-12 < |x| <= 0.62
    FOR(i,to+1) ord_coef[i] = 0;

    enum { ord = DESC_MAX_ORD+1 };
    static num_t scl_coef[ord+1] = {0};
    if (!scl_coef[0]) {
      scl_coef[0] = -1./3;
      FOR(i,1,ord+1) scl_coef[i] = -scl_coef[i-1] * SQR(2*i+1.)/(i*(4*i+6.));
    }

    num_t fact=1;
    num_t tmp_coef[ord+1]; FOR(i,ord+1) tmp_coef[i] = scl_coef[i];
    for (ord_t o = 1; o <= to; o += 2) {
      NUM a0pi = 1, dc1, dc2;
      fact *= o*(o+1.);
      FOR(i,ord+1) {
        tmp_coef[i] *= o > 1 ? -CUB(2.*i+o) / (2*i+o+2) : 1;
        ord_coef [ o          ] += (dc1=a0*SQR(a0pi)*tmp_coef[i]*(  o+1) / fact);
        ord_coef [(o+1)%(to+1)] += (dc2=   SQR(a0pi)*tmp_coef[i]*(2*i+1) / fact);
        if (fabs(dc1)/fabs(ord_coef[o])            < SQR(DBL_EPSILON) &&
            fabs(dc2)/fabs(ord_coef[(o+1)%(to+1)]) < SQR(DBL_EPSILON)) {
          if (DBG_SERIES)
            printf("asinhc: i=%d, o=%d, to=%d, |a0|=%.16e, |dc|=%.16e\n", i, o, to, fabs(a0), MAX(fabs(dc1),fabs(dc2)));
          break;
        } else if (i == ord) {
          if (DBG_SERIES)
            printf("asinhc: i=%d, o=%d, to=%d, |a0|=%.16e, |dc|=%.16e\n", i, o, to, fabs(a0), MAX(fabs(dc1),fabs(dc2)));
          trace(1, "asinhc did not converged after %d iterations", i);
        }
        a0pi *= a0;
      }
    }
    ord_coef[0] = f0;
  } else { // asinhc(x), |x| <= 1e-12
    ord_coef[0] = 1;
    ord_coef[1] = 0;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = -(ord_coef[o-2] * SQR(o-1.)) / (o*(o+1.));
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(erf) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  // erf(z) = 2/sqrt(pi) \int_0^z exp(-t^2) dt
  NUM a0 = a->coef[0];
  NUM f0 = NUMF(erf)(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1], a2 = a0*a0, f1 = M_2_SQRTPI*exp(-a2);
  ord_coef[0] = f0;
  ord_coef[1] = f1;
  for (ord_t o = 2; o <= to; ++o)
    ord_coef[o] = -2*((o-2)*ord_coef[o-2]/(o-1.) + ord_coef[o-1]*a0) / o;

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(erfc) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  FUN(erf)(a,c);
  FUN(axpb)(-1,c,1,c);
  DBGFUN(<-);
}

void
FUN(erfi) (const T *a, T *c) // erfi(z) = -i erf(iz)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
#ifdef MAD_CTPSA_IMPL
  mad_ctpsa_scl(a, I, c);
  mad_ctpsa_erf(c, c);
  mad_ctpsa_scl(c,-I, c);
#else
  ctpsa_t *t = GET_TMPC(c);
  mad_ctpsa_cplx(a , NULL, t);
  mad_ctpsa_scl (t, I, t);
  mad_ctpsa_erf (t, t);
  mad_ctpsa_scl (t,-I, t);
  mad_ctpsa_real(t, c);
  REL_TMPC(t);
#endif
  DBGFUN(<-);
}

void
FUN(erfcx) (const T *a, T *c) // normalized erfc, erfcx(z) = wf(iz)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
#ifdef MAD_CTPSA_IMPL
  mad_ctpsa_scl(a, I, c);
  mad_ctpsa_wf (c, c);
#else
  ctpsa_t *t = GET_TMPC(c);
  mad_ctpsa_cplx(a , NULL, t);
  mad_ctpsa_scl (t, I, t);
  mad_ctpsa_wf  (t, t);
  mad_ctpsa_real(t, c);
  REL_TMPC(t);
#endif
  DBGFUN(<-);
}

void
FUN(wf) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  NUM f0 = NUMF(wf)(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = f0;
  ord_coef[1] = -2*a0*f0 + I*M_2_SQRTPI;

  if (fabs(a0) <= 7 || (cimag(a0) < creal(a0)+0.5 && cimag(a0) < -creal(a0)+0.5)) {
    NUM p[to+1], q[to+1];
    p[0] = 1;
    p[1] = -2*a0;
    q[0] = 0;
    q[1] = 1;

    for (ord_t o = 2; o <= to; ++o) {
      q[o] = 0;
#ifdef MAD_CTPSA_IMPL
      q[o] = -2*(o-1)*q[o-2] - 2*a0*q[o-1];
#endif
      p[o] = -2*(o-1)*p[o-2] - 2*a0*p[o-1];
      ord_coef[o] = (p[o]*f0 + q[o]*I*M_2_SQRTPI)/mad_num_fact(o);
    }

  } else if (fabs(a0) > 7 && fabs(a0) < 100) { // adjusted for double precision
#ifdef MAD_CTPSA_IMPL
    NUM a02 = a0*a0, a04=a02*a02, a06=a04*a02, a08=a04*a04;
    NUM coef[7] = {
      -I*M_SQRTPI/(a04),
      -I*M_SQRTPI/(a06)       *5,
      -I*M_SQRTPI/(a06*a02)/4 *105,
      -I*M_SQRTPI/(a06*a04)/4 *630,
      -I*M_SQRTPI/(a08*a04)/16*17325,
      -I*M_SQRTPI/(a08*a06)/16*135135,
      -I*M_SQRTPI/(a08*a08)/64*4729725 };
#endif

    NUM p[to+1];
    NUM a0r = creal(a0);
    NUM ez2 = exp(-a0r*a0r);
    p[0] = 1;
    p[1] = -2*a0r;
    p[2] = a0r ? 4*a0r*a0r -2 : 0;

    ord_coef[2] = -a0*(-2*a0*f0 + I*M_2_SQRTPI) - f0;

    for (ord_t o = 3; o <= to; ++o) {
      p[o] = -2*(o-1)*p[o-2] - 2*a0r*p[o-1];
      ord_coef[o] = ez2*p[o]/mad_num_fact(o);
#ifdef MAD_CTPSA_IMPL
      FOR(i,7) ord_coef[o] += coef[i];
      FOR(i,7) coef[i] *= -1/a0*(2.*i+o+1)/(o+1);
#endif
    }
  } else {
    for (ord_t o = 1; o <= to; ++o)
      ord_coef[o] = -ord_coef[o-1]/a0;
  }

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(inv_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(inv)(a, CPX(v), c); }

void FUN(invsqrt_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(invsqrt)(a, CPX(v), c); }

void FUN(pown_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(pown)(a, CPX(v), c); }

#endif

// --- end --------------------------------------------------------------------o

/*
Recurrence of Taylor coefficients:
----------------------------------
(f(g(x)))'   = f'(g(x)).g'(x)
(f(x).g(x))' = f'(x)g(x) + f(x)g'(x)

-- sinc(z) -----
[0] sinc(z)     =   sin(z)/z ; cos(z)/z = [c] ; sin(z) = sz ; cos(z) = cz
[1] sinc(z)'    =  cz/z -1sz/z^2                                                     =  [c]-1*[0]/z
[2] sinc(z)''   = -sz/z -2cz/z^2 +2!sz/z^3                                           = -[0]-2*[1]/z
[3] sinc(z)'''  = -cz/z +3sz/z^2 +3!cz/z^3 - 3!sz/z^4                                = -[c]-3*[2]/z
[4] sinc(z)^(4) =  sz/z +4cz/z^2 -12sz/z^3 - 4!cz/z^4 + 4!sz/z^5                     =  [0]-4*[3]/z
[5] sinc(z)^(5) =  cz/z -5sz/z^2 -20cz/z^3 + 60sz/z^4 + 5!cz/z^5 -5!sz/z^6           =  [c]-5*[4]/z
[6] sinc(z)^(6) = -sz/z -6cz/z^2 +30sz/z^3 +120cz/z^4 -360sz/z^5 -6!cz/z^6 +6!sz/z^7 = -[0]-6*[5]/z

-- erf(z) -----
[0] erf(z)
[1] erf(z)'    =     1
[2] erf(z)''   = -2*      z                                   = -2*(0*[0]+[1]*z)
[3] erf(z)'''  = -2*(1 -2*z^2)                                = -2*(1*[1]+[2]*z)
[4] erf(z)^(4) = -2*(-4*z -2*(1-2*z^2)*z)                     = -2*(2*[2]+[3]*z)
[5] erf(z)^(5) = -2*(-6*(1-2*z^2) +4*(3*z-2*z^3)*z)           = -2*(3*[3]+[4]*z)
[6] erf(z)^(6) = -2*(16*(3*z-2*z^3) +4*(3-12*z^2+4*z^4)*z)    = -2*(4*[4]+[5]*z)
                   % *exp(-z^2) *2/sqrt(pi)

demonstration
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

-- erfcx(z) -----
erfcx(z)     = exp(z^2) erfc(z)        =              H_{-1}(z) * 2/sqrt(pi)
erfcx(z)'    = H_{-1}(z)' * 2/sqrt(pi) = -       2    H_{-2}(z) * 2/sqrt(pi)
erfcx(z)''   = H_{-2}(z)' * 2/sqrt(pi) = -     2*4    H_{-3}(z) * 2/sqrt(pi)
erfcx(z)'''  = H_{-3}(z)' * 2/sqrt(pi) = -   2*4*6    H_{-4}(z) * 2/sqrt(pi)
erfcx(z)^(4) = H_{-4}(z)' * 2/sqrt(pi) = - 2*4*6*8    H_{-5}(z) * 2/sqrt(pi)
erfcx(z)^(5) = H_{-5}(z)' * 2/sqrt(pi) = - 2*4*6*8*10 H_{-6}(z) * 2/sqrt(pi)

demonstration
H_n(z)' = 2n H_{n-1}(z)
H_n(z)  = 2z H_{n-1}(z) - 2n H_{n-2}(z)        , n >  0
H_n(z)  = 1/(2n+2) [2z H_{n+1}(z) - H_{n+2}(z)], n < -1

H_n(z)  = 1/(2n+2) [2z H_{n+1}(z) - H_{n+2}(z)], n < -1

H_{ 0}(z) = 1
H_{-1}(z) = erfcx(z) / (2/sqrt(pi))
H_{-2}(z) = -1/2 [2z H_{-1}(z) - H_{ 0}(z)]
H_{-3}(z) = -1/4 [2z H_{-2}(z) - H_{-1}(z)]
H_{-4}(z) = -1/6 [2z H_{-3}(z) - H_{-2}(z)]
H_{-5}(z) = -1/8 [2z H_{-4}(z) - H_{-3}(z)]

n=-(o+1)
H_{n}(z) = H[o+1] = 1/(-2*o) [2z H[o] - H[o-1]] = [0.5*H[o-1] - z H[o]] / o

Plot[HermiteH(-4, z) - (-1/6 (2z HermiteH(-3,z) - HermiteH(-2,z))), {z,-10,10}]

-- Wf(z) -----
Wf(z)     = exp(-z^2) erfc(-iz)
Wf(z)'    =  -2    z^1                       Wf(z) +                     2i/sqrt(pi)
Wf(z)''   =   2  (2z^2 - 1)                  Wf(z) - 2    z^1            2i/sqrt(pi)
Wf(z)'''  =  -4  (2z^3 - 3z)                 Wf(z) + 4  ( z^2- 1)        2i/sqrt(pi)
Wf(z)^(4) =   4  (4z^4 - 12z^2 +  3)         Wf(z) - 4  (2z^3- 5z)       2i/sqrt(pi)
Wf(z)^(5) =  -8  (4z^5 - 20z^3 + 15z)        Wf(z) + 8  (2z^4- 9z^2+ 4)  2i/sqrt(pi)
Wf(z)^(6) =   8  (8z^6 - 60z^4 + 90z^2 - 15) Wf(z) - 8  (4z^5-28z^3+33z) 2i/sqrt(pi)
*/


/* OLD erfcx & wf

void
FUN(erfcx) (const T *a, T *c) // normalized erfc
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  NUM a0 = a->coef[0];
  NUM f0 = NUMF(erfcx)(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1], H[to+2];
  ord_coef[0] = f0;
  H[0] = 1;               // H_0
  H[1] = f0 / M_2_SQRTPI; // H_{-1}
  num_t cc = M_2_SQRTPI;

  for (ord_t o = 1; o <= to; ++o) {
    cc *= -2; // *o => o!
    H[o+1] = (0.5*H[o-1] - a0*H[o]) / o; // H_{-o-1}
    ord_coef[o] = cc*H[o+1]; // / o!
#if 0
#ifdef MAD_CTPSA_IMPL
    printf("o=%d, n=%d, H_n=%+.16e%+.16e, C_o=%+.16e%+.16e\n", o, -(o+1), creal(H[o+1]), cimag(H[o+1]), creal(ord_coef[o]), cimag(ord_coef[o]));
#else
    printf("o=%d, n=%d, H_n=%+.16e, C_o=%+.16e\n", o, -(o+1), H[o+1], ord_coef[o]);
#endif
#endif
  }

// n=-(o+1)
// H_{n}(z) = H[o+1] = 1/(-2*o) [2z H[o] - H[o-1]] = [H[o-1]/2 - zH[o]] / o

  fun_taylor(a,c,to,ord_coef);
  DBGFUN(<-);
}

void
FUN(wf) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");
  NUM a0 = a->coef[0];
  NUM f0 = NUMF(wf)(a0);

  ord_t to = c->mo;
  if (!to || FUN(isval)(a)) { FUN(setval)(c,f0); DBGFUN(<-); return; }

  NUM ord_coef[to+1];

  printf("a0=" FMT "\n", VAL(a0));

  if (fabs(a0) > 10) {
    ord_coef[0] = f0/a0;
    ord_coef[1] = (-2*a0*f0 + I*M_2_SQRTPI)/a0;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = -2*(ord_coef[o-1] + (o-1)*ord_coef[o-2]/a0);
    for (ord_t o = 2; o <= to; ++o) ord_coef[o] /= mad_num_fact(o);
    fun_taylor(a,c,to,ord_coef);
    FUN(scl)(c, a0, c);
  } else {
    ord_coef[0] = f0;
    ord_coef[1] = -2*a0*f0 + I*M_2_SQRTPI;
    for (ord_t o = 2; o <= to; ++o)
      ord_coef[o] = -2*(a0*ord_coef[o-1] + (o-1)*ord_coef[o-2]);
    for (ord_t o = 2; o <= to; ++o) ord_coef[o] /= mad_num_fact(o);
    fun_taylor(a,c,to,ord_coef);
  }
  DBGFUN(<-);
}
*/
