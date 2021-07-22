/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map exponential module implementation
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

#include <math.h>
#include <float.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc (ssz_t sa, const T *ma[sa])
{
  assert(ma);
  for (idx_t i=1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_exppb (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  assert(ma && mb && mc);
  ensure(sa>0 && sb>0, "invalid map sizes (zero or negative sizes)");
  ensure(sb == ma[0]->d->nmv, "incompatibles GTPSA (number of map variables differ)");
  check_same_desc(sa,ma);
  check_same_desc(sb,mb);
  check_same_desc(sa,(const T**)mc);
  ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
print_damap (ssz_t sa, const T *ma[sa], FILE *fp)
{
  char nam[3] = "#1";
  for (ssz_t i=0; i < sa; i++, nam[1]++) FUN(print)(ma[i], nam, 1e-15, 0, fp);
  (void)print_damap;
}

static inline void
exppb1 (ssz_t sa, const T *ma[sa], const T *b, T *c, T *t[3], log_t inv)
{
  FUN(copy)(b, t[0]);
  FUN(copy)(b, c);

  num_t nrm0 = INFINITY;
  num_t nrm_min = 10*DBL_EPSILON*sa;

//fprintf(stderr, "exppb1: X=\n"); print_damap(1, (const T**)&b, stderr);
//fprintf(stderr, "exppb1: H=\n"); print_damap(sa, ma, stderr);

  for (idx_t i=1; i < 100; ++i) { // arbitrary boundary, should be enough!

//fprintf(stderr, "exppb1: i=%d, nrm0=%.10g\n", i, nrm0);
//fprintf(stderr, "exppb1: t[0](b1)=\n"); print_damap(1, (const T**)&t[0], stderr);

    FUN(scl)(t[0], 1.0/i, t[1]);     // t[1] = t[0]/i

//fprintf(stderr, "exppb1: t[1](b2)=\n"); print_damap(1, (const T**)&t[1], stderr);

    // -> c_bra_v_ct
    FUN(clear)(t[0]);

//fprintf(stderr, "bra: n=%d\n", sa);

    for (idx_t j = 0; j < sa; ++j) { // t[0] = h*t[1]
//fprintf(stderr, "bra: j=%d\n", j);

      FUN(deriv)(t[1], t[2], j+1);

//fprintf(stderr, "bra: t[2](s2.d.i)=\n"); print_damap(1, (const T**)&t[2], stderr);

      FUN(mul)(ma[j], t[2], t[3]);
      (inv ? FUN(sub) : FUN(add))(t[0], t[3], t[0]);

//fprintf(stderr, "bra: t[0](s22)=\n"); print_damap(1, (const T**)&t[0], stderr);
    } // <- c_bra_v_ct

//fprintf(stderr, "exppb1: t[0](b1)=\n"); print_damap(1, (const T**)&t[0], stderr);

    FUN(add)(t[0], c, c);           // b3 = t[0]+b3

//fprintf(stderr, "exppb1: c(b3)=\n"); print_damap(1, (const T**)&c, stderr);

    num_t nrm = FUN(nrm)(t[0]);
//fprintf(stderr, "exppb1: nrm0(%d)=%.16e -> nrm(%d)=%.16e\n", i, nrm0, i, nrm);

    if (nrm < nrm_min || nrm >= nrm0) {
//fprintf(stderr, "exppb1: breaking nrm0(%d)=%.16e -> nrm(%d)=%.16e\n", i,nrm0,i,nrm);
      break;
    }
    nrm0 = nrm;
  }
}

// --- public -----------------------------------------------------------------o

void // compute M x = exp(:f(x;0):) x (eq. 32, 33 & 38 and inverse)
FUN(exppb) (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa], int inv)
{
  DBGFUN(->);
  ensure(inv == 1 || inv == -1, "invalid inv value, -1 or 1 expected, got %d", inv);
  check_exppb(sa, ma, sb, mb, mc);

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  for (idx_t ic = 0; ic < sa; ++ic) {
    DBGTPSA(ma[ic]); DBGTPSA(mc[ic]);
    mc_[ic] = FUN(new)(mc[ic], mad_tpsa_same);
  }

  for (idx_t ib = 0; ib < sb; ++ib) {
    DBGTPSA(mb[ib]);
  }

  // temporaries
  T *t[4];
  for (int i = 0; i < 4; ++i) t[i] = FUN(new)(mc[0], mad_tpsa_same);

  for (idx_t i = 0; i < sa; ++i)
    exppb1(sa, ma, mb[i], mc_[i], t, inv == -1);

  // temporaries
  for (int i = 0; i < 4; i++) FUN(del)(t[i]);

  // copy back
  for (idx_t ic = 0; ic < sa; ++ic) {
    FUN(copy)(mc_[ic], mc[ic]);
    FUN(del )(mc_[ic]);
    DBGTPSA(mc[ic]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void // compute G(x;0) = -J grad.f(x;0) (eq. 34),
FUN(vec2fld) (ssz_t sc, const T *a, T *mc[sc]) // cpbbra (wo * -2i)
{
  DBGFUN(->);
  assert(a && mc);
  check_same_desc(sc,(const T**)mc);
  ensure(a->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");

  T *t = FUN(new)(a, mad_tpsa_same);

  for (idx_t ic = 0; ic < sc; ++ic) {
    FUN(setvar)(t, 0, ic+1, 0);
    FUN(poisson)(a, t, mc[ic], 0);
  }

  FUN(del)(t);
  DBGFUN(<-);
}

void // compute f(x;0) = \int_0^x J G(x';0) dx' = x^t J phi G(x;0) (eq. 34, 36 & 37)
FUN(fld2vec) (ssz_t sa, const T *ma[sa], T *c) // cgetpb (wo / -2i)
{
  DBGFUN(->);
  assert(ma && c);
  check_same_desc(sa, ma);
  ensure(ma[0]->d == c->d, "incompatibles GTPSA (descriptors differ)");

  FUN(reset0)(c);

  T *t1 = FUN(new)(c, mad_tpsa_same);
  T *t2 = FUN(new)(c, mad_tpsa_same);

  const idx_t *o2i = c->d->ord2idx;

  for (idx_t ia = 0; ia < sa; ++ia) {
    idx_t iv = ia & 1 ? ia : ia+2;
    FUN(setvar)(t2, 0, iv, 0);
    FUN(mul)(ma[ia], t2, t1);

    ord_t lo = MIN(t1->lo,2); // 2..hi, avoid NaN and Inf
    for (ord_t o = lo; o <= t1->hi; ++o)
      for (idx_t i = o2i[o]; i < o2i[o+1]; ++i)
        t1->coef[i] /= o;

    (ia & 1 ? FUN(add) : FUN(sub))(c, t1, c);
  }

  FUN(del)(t2);
  FUN(del)(t1);
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
