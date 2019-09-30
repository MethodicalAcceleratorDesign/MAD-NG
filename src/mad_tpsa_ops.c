/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA operators module implementation
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

#include "mad_log.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- multiplication helpers -------------------------------------------------o

static inline void
hpoly_diag_mul(const NUM *ca, const NUM *cb, NUM *cc, ssz_t nb,
                 const idx_t l[], const idx_t *idx[])
{
  // asymm: c[2 2] = a[2 0]*b[0 2] + a[0 2]*b[2 0]
  for (idx_t ib = 0; ib < nb; ib++)
    if (cb[ib] || ca[ib])
      for (idx_t ia = idx[0][ib]; ia < idx[1][ib]; ia++) {
        idx_t ic = l[hpoly_idx(ib,ia,nb)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca[ia]*cb[ib] + (ia == ib ? 0 : ca[ib]*cb[ia]);
      }
}

static inline void
hpoly_sym_mul(const NUM *ca1, const NUM *cb1, const NUM *ca2, const NUM *cb2,
              NUM *cc, ssz_t na, ssz_t nb, const idx_t l[], const idx_t *idx[])
{
  // na > nb so longer loop is inside
  for (idx_t ib=0; ib < nb; ib++)
    if (cb1[ib] || ca2[ib])
      for (idx_t ia=idx[0][ib]; ia < idx[1][ib]; ia++) {
        idx_t ic = l[hpoly_idx(ib, ia, na)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca1[ia]*cb1[ib] + ca2[ib]*cb2[ia];
      }
}

static inline void
hpoly_asym_mul(const NUM *ca, const NUM *cb, NUM *cc, ssz_t na, ssz_t nb,
               const idx_t l[], const idx_t *idx[])
{
  // oa > ob so longer loop is inside
  for (idx_t ib=0; ib < nb; ib++)
    if (cb[ib])
      for (idx_t ia=idx[0][ib]; ia < idx[1][ib]; ia++) {
        idx_t ic = l[hpoly_idx(ib,ia,na)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca[ia]*cb[ib];
      }
}

static inline void
hpoly_mul(const T *a, const T *b, T *c, const ord_t *ocs, bit_t *cnz, log_t in_parallel)
{
  const D *d = c->d;
  const NUM *ca = a->coef, *cb = b->coef;
  NUM *cc = c->coef;
  idx_t *pi = d->ord2idx, hod = d->mo/2;
  bit_t nza = a->nz, nzb = b->nz;
  ord_t lo = MAX(c->lo,3);

//  printf("%s(%d,%d->%d)\n", __func__, in_parallel, ocs[0], lo);
  for (ord_t i = 0; ocs[i] >= lo; ++i) {
    ord_t oc = ocs[i];
    idx_t idx0 = 0, idx1 = 2;

    if (in_parallel && ocs[i] >= c->hi) {
      oc = c->hi;
      if (ocs[i] == c->hi) idx1 = 1;
      else                 idx0 = 1;
    }

    for (ord_t j=1; j <= (oc-1)/2; ++j) {
      ord_t oa = oc-j, ob = j;            // oa > ob >= 1
      ssz_t na = pi[oa+1] - pi[oa];
      ssz_t nb = pi[ob+1] - pi[ob];
      const idx_t *lc = d->L[oa*hod + ob];
      const idx_t *idx[2] = { d->L_idx[oa*hod + ob][idx0],
                              d->L_idx[oa*hod + ob][idx1] };
       assert(lc); assert(idx[0] && idx[1]);

      if (mad_bit_get(nza & nzb,oa) && mad_bit_get(nza & nzb,ob)) {
//        printf("hpoly__sym_mul (%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_sym_mul(ca+pi[oa],cb+pi[ob], ca+pi[ob],cb+pi[oa], cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_get(nza,oa) && mad_bit_get(nzb,ob)) {
//        printf("hpoly_asym_mul1(%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_asym_mul(ca+pi[oa],cb+pi[ob],cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_get(nza,ob) && mad_bit_get(nzb,oa)) {
//        printf("hpoly_asym_mul2(%d) %2d+%2d=%2d\n", ocs[0], ob,oa,oc);
        hpoly_asym_mul(cb+pi[oa],ca+pi[ob],cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
    }
    // even oc, diagonal case
    if (!(oc & 1) && mad_bit_get(nza & nzb,oc/2)) {
      ord_t hoc = oc/2;
      ssz_t nb = pi[hoc+1]-pi[hoc];
      const idx_t *lc = d->L[hoc*hod + hoc];
      const idx_t *idx[2] = { d->L_idx[hoc*hod + hoc][idx0],
                              d->L_idx[hoc*hod + hoc][idx1] };
      assert(lc); assert(idx[0] && idx[1]);
//      printf("hpoly_diag_mul (%d) %2d+%2d=%2d\n", ocs[0], hoc,hoc,oc);
      hpoly_diag_mul(ca+pi[hoc],cb+pi[hoc],cc, nb, lc, idx);
      *cnz = mad_bit_set(*cnz,oc);
    }
  }
}

#ifdef _OPENMP
static inline void
hpoly_mul_par(const T *a, const T *b, T *c) // parallel version
{
  const D *d = c->d;
  bit_t c_nzs[d->nth];

//  printf("%s(%d)\n", __func__, c->d->nth);
//  FUN(debug)(a,"ain",0);
//  FUN(debug)(b,"bin",0);
//  FUN(debug)(c,"cin",0);

#pragma omp parallel for
  for (int t = 0; t < d->nth; ++t) {
    c_nzs[t] = c->nz;
    ord_t i = 0; while (d->ocs[1+t][i] > c->hi+1) ++i;
    hpoly_mul(a, b, c, &d->ocs[1+t][i], &c_nzs[t], TRUE);
  }

  for (int t = 0; t < d->nth; ++t)
    c->nz |= c_nzs[t];

//  FUN(debug)(c,"cout",0);
}
#endif

static inline void
hpoly_mul_ser(const T *a, const T *b, T *c) // serial version
{
//  printf("%s(1)\n", __func__);
//  FUN(debug)(a,"ain",0);
//  FUN(debug)(b,"bin",0);
//  FUN(debug)(c,"cin",0);

  hpoly_mul(a, b, c, &c->d->ocs[0][c->d->mo-c->hi], &c->nz, FALSE);

//  FUN(debug)(c,"cout",0);
}

// --- derivative helpers -----------------------------------------------------o

static inline int
der_coef(idx_t ia, idx_t di, ord_t der_ord, const D* d)
{
  if (der_ord == 1) // der against var
    return d->To[ia][di-1];

  const ord_t *srcm = d->To[ia], *derm = d->To[di];
  if (mad_mono_gt(d->nv,derm,srcm))
    return 0;

  int c = 1;
  for (int v = 0; v < d->nv; ++v)
    for (int o = 0; o < derm[v]; ++o)
        c *= srcm[v] - o;
  return c;
}

static inline void
hpoly_der_lt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
  const ord_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *pi = d->ord2idx;
  idx_t nc = pi[oc+1] - pi[oc], cols = pi[ord+1] - pi[ord];
  // idx = idx - pi[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(ic,idx-pi[ord],cols)];
    if (ia >= 0 && ca[ia]) {
      assert(pi[oc+ord] <= ia && ia < pi[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    }
  }
}

static inline void
hpoly_der_eq(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
  const ord_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *pi = d->ord2idx;
  idx_t nc = pi[ord+1] - pi[ord];
  idx_t idx_shifted = idx - pi[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(MAX(ic,idx_shifted),MIN(ic,idx_shifted),nc)];
    if (ia >= 0 && ca[ia]) {
      assert(pi[oc+ord] <= ia && ia < pi[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    }
  }
}

static inline void
hpoly_der_gt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
  const ord_t ho = d->mo/2;
  const idx_t *lc = d->L[oc*ho + ord], *pi = d->ord2idx;
  idx_t nc = pi[oc+1] - pi[oc];
  idx_t idx_shifted = idx - pi[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(idx_shifted,ic,nc)];
    if (ia >= 0 && ca[ia]) {
      assert(pi[oc+ord] <= ia && ia < pi[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    }
  }
}

static inline void
hpoly_der(const T *a, idx_t idx, ord_t ord, T *c)
{
  const D *d = c->d;
  idx_t *pi = d->ord2idx;
  const NUM *ca = a->coef;
        NUM *cc;

  c->hi = MIN3(c->mo, d->to, a->hi-ord);  // initial guess, readjust based on nz
  for (ord_t oc = 1; oc <= c->hi; ++oc)
    if (mad_bit_get(a->nz,oc+ord)) {
      cc = c->coef + pi[oc];
           if (oc <  ord) hpoly_der_lt(ca,cc,idx,oc,ord,&c->nz,d);
      else if (oc == ord) hpoly_der_eq(ca,cc,idx,oc,ord,&c->nz,d);
      else                hpoly_der_gt(ca,cc,idx,oc,ord,&c->nz,d);
    }
  ord_t n = mad_bit_lowest(c->nz);
  c->lo = MIN(n,c->mo);
  c->hi = mad_bit_highest(c->nz);
  if (c->lo) c->coef[0] = 0;
}

// --- binary ops -------------------------------------------------------------o

// TPSA_LINOP(+, +, 0) => cc[i] = +ca[i] + cb[i], with i from (lo+0) to hi

#define TPSA_LINOP(OPA, OPB, ORD) \
do { \
  idx_t *pi = c->d->ord2idx; \
  idx_t start_a = pi[ORD?MAX(a->lo,ORD):a->lo], end_a = pi[MIN(a->hi,c_hi)+1]; \
  idx_t start_b = pi[ORD?MAX(b->lo,ORD):b->lo], end_b = pi[MIN(b->hi,c_hi)+1]; \
  idx_t i = start_a; \
  for (; i < MIN(end_a,start_b); ++i) c->coef[i] = OPA a->coef[i]; \
  for (; i <           start_b ; ++i) c->coef[i] = 0; \
  for (; i < MIN(end_a,end_b)  ; ++i) c->coef[i] = OPA a->coef[i] OPB b->coef[i]; \
  for (; i <     end_a         ; ++i) c->coef[i] = OPA a->coef[i]; \
  for (; i <           end_b   ; ++i) c->coef[i] =                OPB b->coef[i]; \
} while(0)

#define TPSA_CLEAR(ORD) \
do { \
  idx_t *pi = c->d->ord2idx; \
  for (idx_t i = pi[ORD]; i < pi[c_hi+1]; ++i) c->coef[i] = 0; \
} while(0)

void
FUN(scl) (const T *a, NUM v, T *c)
{
  assert(a && c);
  const D *d = a->d;
  ensure(d == c->d, "incompatibles GTPSA (descriptors differ)");

  if (!v || a->hi == 0) { FUN(scalar)(c, v*a->coef[0], 0, 0); return; }

  FUN(copy0)(a,c);

  idx_t *pi = d->ord2idx;
  if (v == 1)
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = a->coef[i];
  else if (v == -1)
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = -a->coef[i];
  else
    for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
      c->coef[i] = v * a->coef[i];
  CHECK_VALIDITY(c);
}

void
FUN(acc) (const T *a, NUM v, T *c)
{
  assert(a && c);
  const D *d = c->d;
  ensure(a->d == d, "incompatibles GTPSA (descriptors differ)");

  if (!v) return;

  ord_t new_lo = MIN (a->lo,c->lo);
  ord_t new_hi = MIN3(a->hi,c->mo,d->to);
  idx_t *pi = d->ord2idx;

  for (idx_t i = pi[new_lo ]; i < pi[c->lo   ]; ++i) c->coef[i]  = 0;
  for (idx_t i = pi[c->hi+1]; i < pi[new_hi+1]; ++i) c->coef[i]  = 0;
  for (idx_t i = pi[a->lo  ]; i < pi[new_hi+1]; ++i) c->coef[i] += v * a->coef[i];

  c->lo = new_lo;
  c->hi = MAX(new_hi, c->hi);
  c->nz = mad_bit_hcut(c->nz|a->nz,c->hi);
  if (c->lo) c->coef[0] = 0;

  CHECK_VALIDITY(c);
}

void
FUN(add) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  const D *d = a->d;
  ensure(d == b->d && d == c->d, "incompatibles GTPSA (descriptors differ)");

  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, d->to);
  TPSA_LINOP(, +, 0);  // c->coef[i] = a->coef[i] + b->coef[i];

  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_hcut(a->nz|b->nz, c_hi);
  if (c->lo) c->coef[0] = 0;

  CHECK_VALIDITY(c);
}

void
FUN(sub) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  const D *d = a->d;
  ensure(d == b->d && d == c->d, "incompatibles GTPSA (descriptors differ)");

  const T* t = 0;
  if (a->lo > b->lo) SWAP(a,b,t);

  ord_t   hi = MAX (a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, d->to);
  if (t) TPSA_LINOP(-, +, 0); // c->coef[i] = - a->coef[i] + b->coef[i];
  else   TPSA_LINOP( , -, 0); // c->coef[i] =   a->coef[i] - b->coef[i];

  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_hcut(a->nz|b->nz, c_hi);
  if (c->lo) c->coef[0] = 0;

  CHECK_VALIDITY(c);
}

void
FUN(mul) (const T *a, const T *b, T *r)
{
  assert(a && b && r);
  const D *d = a->d;
  ensure(d == b->d && d == r->d, "incompatibles GTPSA (descriptors differ)");

  T *c = (a == r || b == r) ? GET_TMPX(r) : r;

  c->lo = a->lo + b->lo;
  c->hi = MIN3(a->hi + b->hi, c->mo, d->to);
  // empty
  if (c->lo > c->hi) { FUN(reset0)(c); goto ret; }
  // a is the left-most one
  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }

  // order 0
  NUM a0 = a->coef[0], b0 = b->coef[0];
  c->coef[0] = a0 * b0;
  c->nz = c->coef[0] != 0;
  if (c->hi == 0) {
    if (!c->coef[0]) c->lo = c->mo; // reset
    goto ret;
  }

  // order 1
  if (c->lo <= 1) {
    idx_t max_ord1 = d->ord2idx[2];
    if (mad_bit_get(a->nz & b->nz,1) && a0 && b0) {
      for (idx_t i = 1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i] + b0*a->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
    else if (mad_bit_get(a->nz,1) && b0) {
      for (idx_t i = 1; i < max_ord1; ++i) c->coef[i] = b0*a->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
    else if (mad_bit_get(b->nz,1) && a0) {
      for (idx_t i = 1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
  }

  // order 2+
  if (c->hi > 1) {
    ord_t c_hi = c->hi, ab_hi = MAX(a->hi,b->hi);
    if (a0 && b0) TPSA_LINOP(b0*, +a0*, 2);
    else if (!b0) TPSA_LINOP( 0*, +a0*, 2);
    else          TPSA_LINOP(b0*, + 0*, 2);
    TPSA_CLEAR(ab_hi+1);

    if (a0) { c->nz |= mad_bit_lcut(b->nz,2); }
    if (b0) { c->nz |= mad_bit_lcut(a->nz,2); }
              c->nz  = mad_bit_hcut(c->nz,c_hi);

    if (mad_bit_get(a->nz & b->nz,1)) {
      idx_t *pi = d->ord2idx, hod = d->mo/2;
      const idx_t *lc = d->L[hod+1];
      const idx_t *idx[2] = { d->L_idx[hod+1][0], d->L_idx[hod+1][2] };
      assert(lc);
      hpoly_diag_mul(a->coef+pi[1], b->coef+pi[1], c->coef, pi[2]-pi[1], lc, idx);
      c->nz = mad_bit_set(c->nz,2);
    }

    // order 3+
    if (c->hi > 2) {
#ifdef _OPENMP
#ifdef MAD_TPSA_MUL_PAR
      //if (c->hi >= 12)
      if (d->ord2idx[c->hi+1] >= 1000000) // parallel mul is always slower...
        hpoly_mul_par(a,b,c);
      else
#endif
#endif
        hpoly_mul_ser(a,b,c);
    }
  }

ret:
  assert(a != c && b != c);
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  CHECK_VALIDITY(c);
}

void
FUN(div) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  NUM b0 = b->coef[0];
  ensure(b0 != 0, "invalid domain");

  if (b->hi == 0) { FUN(scl) (a,1/b0,c); return; }

  T *t = (a == c || b == c) ? GET_TMPX(c) : c;
  FUN(inv)(b,1,t);
  FUN(mul)(a,t,c);
  if (t != c) REL_TMPX(t);

  CHECK_VALIDITY(c);
}

void
FUN(powi) (const T *a, int n, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  int inv = 0;

  if (n < 0) { n = -n; inv = 1; }

  T *t1 = GET_TMPX(c);

  switch (n) {
    case 0: FUN(scalar) (c, 1, 0, 0); break; // ok: no copy
    case 1: FUN(copy  ) (a, c);       break; // ok: 1 copy
    case 2: FUN(mul   ) (a,a, c);     break; // ok: 1 copy if a==c
    case 3: FUN(mul   ) (a,a, t1); FUN(mul)(t1,a,  c); break; // ok: 1 copy if a==c
    case 4: FUN(mul   ) (a,a, t1); FUN(mul)(t1,t1, c); break; // ok: no copy
    default: {
      T *t2 = GET_TMPX(c), *t;
      FUN(copy  )(a, t1);
      FUN(scalar)(c, 1, 0, 0);
      for (;;) {
        if (n  & 1)   FUN(mul)(c ,t1, c ); // ok: 1 copy
        if (n /= 2) { FUN(mul)(t1,t1, t2); SWAP(t1,t2,t); } // ok: no copy
        else break;
      }
      REL_TMPX(t2);
    }
  }
  REL_TMPX(t1);

  if (inv) FUN(inv)(c,1,c);
  CHECK_VALIDITY(c);
}

log_t
FUN(equ) (const T *a, const T *b, num_t eps_)
{
  assert(a && b);
  ensure(a->d == b->d, "incompatibles GTPSA (descriptors differ)");

  if (eps_ < 0) eps_ = 1e-16;

  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }

  idx_t *pi = a->d->ord2idx;
  idx_t start_a = pi[a->lo], end_a = pi[a->hi+1];
  idx_t start_b = pi[b->lo], end_b = pi[b->hi+1];
  idx_t i = start_a;

  for (; i < MIN(end_a,start_b); ++i) if (fabs(a->coef[i]) > eps_)              return FALSE;
  i = start_b;
  for (; i < MIN(end_a,end_b)  ; ++i) if (fabs(a->coef[i] - b->coef[i]) > eps_) return FALSE;
  for (; i <     end_a         ; ++i) if (fabs(a->coef[i]) > eps_)              return FALSE;
  for (; i <           end_b   ; ++i) if (fabs(b->coef[i]) > eps_)              return FALSE;

  return TRUE;
}

// --- function ---------------------------------------------------------------o

void
FUN(abs) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  c->lo = a->lo;
  c->hi = MIN3(a->hi, c->mo, c->d->to);
  c->nz = mad_bit_hcut(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = fabs(a->coef[i]);

  CHECK_VALIDITY(c);
}

#ifdef MAD_CTPSA_IMPL

void
FUN(arg) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  c->lo = a->lo;
  c->hi = MIN3(a->hi, c->mo, c->d->to);
  c->nz = mad_bit_hcut(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = carg(a->coef[i]);
}

void
FUN(conj) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  c->lo = a->lo;
  c->hi = MIN3(a->hi, c->mo, c->d->to);
  c->nz = mad_bit_hcut(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (idx_t i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = conj(a->coef[i]);

  CHECK_VALIDITY(c);
}

#endif

NUM
FUN(nrm1) (const T *a, const T *b_)
{
  assert(a);
  NUM norm = 0.0;
  idx_t *pi = a->d->ord2idx;
  if (b_) {
    ensure(a->d == b_->d, "incompatibles GTPSA (descriptors differ)");
    if (a->lo > b_->lo) { const T *t; SWAP(a,b_,t); }

    idx_t start_a = pi[a ->lo], end_a = pi[MIN(a ->hi,a ->d->to)+1],
          start_b = pi[b_->lo], end_b = pi[MIN(b_->hi,b_->d->to)+1];
    idx_t i;
    for (i = start_a; i < MIN(end_a,start_b); ++i) norm += fabs(a->coef[i]);
    for (i = start_b; i < MIN(end_a,end_b)  ; ++i) norm += fabs(a->coef[i] - b_->coef[i]);
    for (           ; i <     end_a         ; ++i) norm += fabs(a->coef[i]);
    for (           ; i <     end_b         ; ++i) norm += fabs(b_->coef[i]);
  }
  else {
    ord_t hi = MIN(a->hi, a->d->to);
    for (ord_t o = a->lo; o <= hi; ++o) {
      if (!mad_bit_get(a->nz,o)) continue;
      for (idx_t i = pi[o]; i < pi[o+1]; ++i)
        norm += fabs(a->coef[i]);
    }
  }
  return norm;
}

NUM
FUN(nrm2) (const T *a, const T *b_)
{
  assert(a);
  NUM norm = 0.0;
  idx_t *pi = a->d->ord2idx;
  if (b_) {
    ensure(a->d == b_->d, "incompatibles GTPSA (descriptors differ)");
    if (a->lo > b_->lo) { const T* t; SWAP(a,b_,t); }

    idx_t start_a = pi[a ->lo], end_a = pi[MIN(a ->hi,a ->d->to)+1],
          start_b = pi[b_->lo], end_b = pi[MIN(b_->hi,b_->d->to)+1];
    idx_t i;
    for (i = start_a; i < MIN(end_a,start_b); ++i) norm +=  a->coef[i]              *  a->coef[i];
    for (i = start_b; i < MIN(end_a,end_b)  ; ++i) norm += (a->coef[i]-b_->coef[i]) * (a->coef[i]-b_->coef[i]);
    for (           ; i <     end_a         ; ++i) norm +=  a->coef[i]              *  a->coef[i];
    for (           ; i <     end_b         ; ++i) norm +=  b_->coef[i]             *  b_->coef[i];
  }
  else {
    ord_t hi = MIN(a->hi, a->d->to);
    for (ord_t o = a->lo; o <= hi; ++o)
      if (mad_bit_get(a->nz,o)) {
        for (idx_t i = pi[o]; i < pi[o+1]; ++i)
          norm += a->coef[i] * a->coef[i];
      }
  }
  return sqrt(norm);
}

void
FUN(deriv) (const T *a, T *c, int iv)
{
  assert(a && c && a != c); // TODO: aliasing
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ensure(iv >= a->d->ord2idx[1] && iv < a->d->ord2idx[2], "invalid domain");
  // TODO: ensure map_order[iv] > 0

  if (!a->hi) { FUN(reset0)(c); return; }
  FUN(scalar)(c,FUN(geti)(a,iv), 0, 0);  // TODO: what if alpha[iv] == 0 ?

  const D *d = c->d;
  c->lo = a->lo ? a->lo-1 : 0;  // initial guess, readjusted after computation
  c->hi = MIN3(a->hi-1, c->mo, d->to);

  idx_t *pi = d->ord2idx;
  const NUM *ca = a->coef;

  ord_t der_ord = 1, oc = 1;
  if (mad_bit_get(a->nz,oc+1))
      hpoly_der_eq(ca,c->coef+pi[oc],iv,oc,der_ord,&c->nz,d);
  for (oc = 2; oc <= c->hi; ++oc)
    if (mad_bit_get(a->nz,oc+1))
      hpoly_der_gt(ca,c->coef+pi[oc],iv,oc,der_ord,&c->nz,d);
  ord_t n = mad_bit_lowest(c->nz);
  c->lo = MIN(n,c->mo);
  c->hi = mad_bit_highest(c->nz);

  CHECK_VALIDITY(c);
}

void
FUN(derivm) (const T *a, T *c, ssz_t n, const ord_t mono[n])
{
  assert(a && c && a != c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");
  ensure(mad_desc_mono_isvalid_m(a->d,n,mono), "invalid monomial");

  ord_t der_ord = mad_mono_ord(n,mono);
  ensure(der_ord > 0, "invalid derivative order");
  idx_t idx = mad_desc_get_idx_m(a->d,n,mono);
  if (idx < a->d->ord2idx[2]) {  // fallback on simple version
    FUN(deriv)(a,c,idx);
    return;
  }

  // ord 0 & setup
  FUN(scalar)(c,FUN(geti)(a,idx) * der_coef(idx,idx,der_ord,a->d), 0, 0);
  if (a->hi <= der_ord)
    return;

  // ords 1..a->hi - 1
  hpoly_der(a,idx,der_ord,c);

  CHECK_VALIDITY(c);
}

void
FUN(poisson) (const T *a, const T *b, T *c, int nv)
{
  // C = [A,B]
  assert(a && b && c);
  ensure(a->d == b->d && b->d == c->d, "incompatibles GTPSA (descriptors differ)");

  nv = nv>0 ? nv/2 : a->d->nv/2;

  T *is[4];
  for (int i = 0; i < 4; ++i)
    is[i] = FUN(new)(a, a->d->to);

  for (int i = 1; i <= nv; ++i) {
    FUN(deriv)(a, is[0], 2*i - 1); // res = res + da/dq_i * db/dp_i
    FUN(deriv)(b, is[1], 2*i    );
    FUN(mul)  (is[0],is[1],is[2]);
    FUN(add)  (is[3],is[2],is[0]);
    FUN(copy) (is[0],is[3]);

    FUN(deriv)(a, is[0], 2*i    ); // res = res - da/dp_i * db/dq_i
    FUN(deriv)(b, is[1], 2*i - 1);
    FUN(mul)  (is[0],is[1],is[2]);
    FUN(sub)  (is[3],is[2],is[0]);
    FUN(copy) (is[0],is[3]);
  }

  FUN(copy)(is[3], c);
  for (int i = 0; i < 4; ++i) FUN(del)(is[i]);
}

// --- high level functions ---------------------------------------------------o

void
FUN(axpb) (NUM a, const T *x, NUM b, T *r)
{
  assert(x && r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");
  FUN(scl)(x,a,r);
  if (b) FUN(set0)(r,1,b);
}

void
FUN(axpbypc) (NUM c1, const T *a, NUM c2, const T *b, NUM c3, T *c)
{            //    a           x       b           y      c      r
  assert(a && b && c);
  ensure(a->d == b->d && b->d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->to);

  if (!c_hi) { FUN(scalar)(c, c1*a->coef[0] + c2*b->coef[0] + c3, 0, 0); return; }

  if (a->lo > b->lo)  {
    const T* t; SWAP(a ,b ,t);
    NUM n;      SWAP(c1,c2,n);
  }

  if (c1 == 1) {
         if (c2 ==  1) TPSA_LINOP( , +   , 0);
    else if (c2 == -1) TPSA_LINOP( , -   , 0);
    else               TPSA_LINOP( , +c2*, 0);
  } else
  if (c1 == -1) {
         if (c2 ==  1) TPSA_LINOP( -, +   , 0);
    else if (c2 == -1) TPSA_LINOP( -, -   , 0);
    else               TPSA_LINOP( -, +c2*, 0);
  } else {
         if (c2 ==  1) TPSA_LINOP( c1*, +   , 0);
    else if (c2 == -1) TPSA_LINOP( c1*, -   , 0);
    else               TPSA_LINOP( c1*, +c2*, 0);
  }

  c->lo = a->lo;    // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_hcut(a->nz|b->nz,c->hi);
  if (c->lo) c->coef[0] = 0;

  if (c3) FUN(set0)(c,1,c3);

  CHECK_VALIDITY(c);
}

void
FUN(axypb) (NUM a, const T *x, const T *y, NUM b, T *r)
{
  assert(x && y && r);
  ensure(x->d == y->d && y->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(mul)(x,y,r);
  if (a != 1 || b != 0) FUN(axpb)(a,r,b,r);
}

void
FUN(axypbzpc) (NUM a, const T *x, const T *y, NUM b, const T *z, NUM c, T *r)
{
  assert(x && y && z && r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = z == r ? GET_TMPX(r) : r;
  FUN(mul)(x,y,t);
  FUN(axpbypc)(a,t,b,z,c,r);
  if (t != r) REL_TMPX(t);
}

void
FUN(axypbvwpc) (NUM a, const T *x, const T *y,
                NUM b, const T *v, const T *w, NUM c, T *r)
{
  assert(x && y && v && w && r);
  ensure(x->d == y->d && y->d == v->d && v->d == w->d && w->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(mul)(x,y,t);
  FUN(mul)(v,w,r);
  FUN(axpbypc)(a,t,b,r,c,r);
  REL_TMPX(t);
}

void
FUN(ax2pby2pcz2) (NUM a, const T *x, NUM b, const T *y, NUM c, const T *z, T *r)
{
  assert(x && y && z && r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = z == r ? GET_TMPX(r) : r;
  FUN(axypbvwpc)(a,x,x, b,y,y, 0, t);
  FUN(axypbzpc)(c,z,z, 1,t, 0, r);
  if (t != r) REL_TMPX(t);
}

void
FUN(axpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");

  T *t = x == r ? GET_TMPX(r) : r;
  FUN(axypb)(c,x,x,b,t);
  FUN(sqrt)(t,t);
  FUN(axpbypc)(a,x,1,t,0,r);
  if (t != r) REL_TMPX(t);
}

void
FUN(logaxpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(axpsqrtbpcx2)(x, a, b, c, r);
  FUN(log)(r, r);
}

void
FUN(logxdy) (const T *x, const T *y, T *r)
{
  assert(x && y && r);
  ensure(x->d == y->d && y->d == r->d, "incompatibles GTPSA (descriptors differ)");

  NUM x0 = FUN(get0)(x), y0 = FUN(get0)(y);
  if (fabs(x0) > fabs(y0)) { // TODO: improve stability around 1i
    FUN(div)(x, y, r);
    FUN(log)(r, r);
  } else {
    FUN(div)(y, x, r);
    FUN(log)(r, r);
    FUN(scl)(r, -1, r);
  }
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(nrm1_r) (const T *t, const T *t2_, cnum_t *r)
{ *r = FUN(nrm1)(t, t2_); }

void FUN(nrm2_r) (const T *t, const T *t2_, cnum_t *r)
{ *r = FUN(nrm2)(t, t2_); }

void FUN(acc_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(acc)(a, CNUM(v), c); }

void FUN(scl_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(scl)(a, CNUM(v), c); }

void FUN(axpb_r) (num_t a_re, num_t a_im, const T *x,
                  num_t b_re, num_t b_im, T *r)
{ FUN(axpb)(CNUM(a), x, CNUM(b), r); }

void FUN(axpbypc_r) (num_t a_re, num_t a_im, const T *x,
                     num_t b_re, num_t b_im, const T *y,
                     num_t c_re, num_t c_im, T *r)
{ FUN(axpbypc)(CNUM(a), x, CNUM(b), y, CNUM(c), r); }

void FUN(axypb_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                   num_t b_re, num_t b_im, T *r)
{ FUN(axypb)(CNUM(a), x, y, CNUM(b), r); }

void FUN(axypbzpc_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                      num_t b_re, num_t b_im, const T *z,
                      num_t c_re, num_t c_im, T *r)
{ FUN(axypbzpc)(CNUM(a), x, y, CNUM(b), z, CNUM(c), r); }

void FUN(axypbvwpc_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                       num_t b_re, num_t b_im, const T *v, const T *w,
                       num_t c_re, num_t c_im, T *r)
{ FUN(axypbvwpc)(CNUM(a), x, y, CNUM(b), v, w, CNUM(c), r); }

void FUN(ax2pby2pcz2_r)(num_t a_re, num_t a_im, const T *x,
                        num_t b_re, num_t b_im, const T *y,
                        num_t c_re, num_t c_im, const T *z, T *r)
{ FUN(ax2pby2pcz2)(CNUM(a), x, CNUM(b), y, CNUM(c), z, r); }

void FUN(axpsqrtbpcx2_r)(const T *x, num_t a_re, num_t a_im,
                                     num_t b_re, num_t b_im,
                                     num_t c_re, num_t c_im, T *r)
{ FUN(axpsqrtbpcx2)(x, CNUM(a), CNUM(b), CNUM(c), r); }

void FUN(logaxpsqrtbpcx2_r) (const T *x, num_t a_re, num_t a_im,
                                         num_t b_re, num_t b_im,
                                         num_t c_re, num_t c_im, T *r)
{ FUN(logaxpsqrtbpcx2)(x, CNUM(a), CNUM(b), CNUM(c), r); }

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
