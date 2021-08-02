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
#include "mad_cst.h"
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
  const idx_t *o2i = d->ord2idx;
  const NUM *ca = a->coef, *cb = b->coef;
  NUM *cc = c->coef;
  idx_t hod = d->mo/2;
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
      ssz_t na = o2i[oa+1] - o2i[oa];
      ssz_t nb = o2i[ob+1] - o2i[ob];
      const idx_t *lc = d->L[oa*hod + ob];
      const idx_t *idx[2] = { d->L_idx[oa*hod + ob][idx0],
                              d->L_idx[oa*hod + ob][idx1] };
       assert(lc); assert(idx[0] && idx[1]);

      if (mad_bit_tst(nza & nzb,oa) && mad_bit_tst(nza & nzb,ob)) {
//        printf("hpoly__sym_mul (%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_sym_mul(ca+o2i[oa],cb+o2i[ob],ca+o2i[ob],cb+o2i[oa],cc,na,nb,lc,idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_tst(nza,oa) && mad_bit_tst(nzb,ob)) {
//        printf("hpoly_asym_mul1(%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_asym_mul(ca+o2i[oa],cb+o2i[ob],cc,na,nb,lc,idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_tst(nza,ob) && mad_bit_tst(nzb,oa)) {
//        printf("hpoly_asym_mul2(%d) %2d+%2d=%2d\n", ocs[0], ob,oa,oc);
        hpoly_asym_mul(cb+o2i[oa],ca+o2i[ob],cc,na,nb,lc,idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
    }
    // even oc, diagonal case
    if (!(oc & 1) && mad_bit_tst(nza & nzb,oc/2)) {
      ord_t hoc = oc/2;
      ssz_t nb = o2i[hoc+1] - o2i[hoc];
      const idx_t *lc = d->L[hoc*hod + hoc];
      const idx_t *idx[2] = { d->L_idx[hoc*hod + hoc][idx0],
                              d->L_idx[hoc*hod + hoc][idx1] };
      assert(lc); assert(idx[0] && idx[1]);
//      printf("hpoly_diag_mul (%d) %2d+%2d=%2d\n", ocs[0], hoc,hoc,oc);
      hpoly_diag_mul(ca+o2i[hoc],cb+o2i[hoc],cc,nb,lc,idx);
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
//  FUN(debug)(a,"ain",0,0);
//  FUN(debug)(b,"bin",0,0);
//  FUN(debug)(c,"cin",0,0);

#pragma omp parallel for
  for (int t = 0; t < d->nth; ++t) {
    c_nzs[t] = c->nz;
    ord_t i = 0; while (d->ocs[1+t][i] > c->hi+1) ++i;
    hpoly_mul(a, b, c, &d->ocs[1+t][i], &c_nzs[t], TRUE);
  }

  for (int t = 0; t < d->nth; ++t)
    c->nz |= c_nzs[t];

//  FUN(debug)(c,"cout",0,0);
}
#endif

static inline void
hpoly_mul_ser(const T *a, const T *b, T *c) // serial version
{
//  printf("%s(1)\n", __func__);
//  FUN(debug)(a,"ain",0,0);
//  FUN(debug)(b,"bin",0,0);
//  FUN(debug)(c,"cin",0,0);

  hpoly_mul(a, b, c, &c->d->ocs[0][c->d->mo-c->hi], &c->nz, FALSE);

//  FUN(debug)(c,"cout",0,0);
}

// --- derivative helpers -----------------------------------------------------o

static inline num_t
der_coef(idx_t ia, idx_t idx, ord_t ord, const D* d)
{
  if (ord == 1) { // der against var
//  fprintf(stderr, "der_coef[%d][%d]=%d\n", ia, idx-1, d->To[ia][idx-1]);
    return d->To[ia][idx-1];
  }
  const ord_t *srcm = d->To[ia], *derm = d->To[idx];
  if (mad_mono_gt(d->nv,derm,srcm))
    return 0;

  num_t c = 1;
  for (int v = 0; v < d->nv; ++v)
    for (ord_t o = 0; o < derm[v]; ++o)
      c *= srcm[v] - o;
  return c;
}

static inline void // oc < ord
hpoly_der_lt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
//fprintf(stderr, "hpoly_der_lt: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *o2i = d->ord2idx;
  idx_t nc = o2i[oc+1] - o2i[oc], cols = o2i[ord+1] - o2i[ord];
  idx_t idx_shifted = idx - o2i[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(ic,idx_shifted,cols)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    } else cc[ic] = 0;
  }
}

static inline void // oc == ord
hpoly_der_eq(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
//fprintf(stderr, "hpoly_der_eq: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *o2i = d->ord2idx;
  idx_t nc = o2i[ord+1] - o2i[ord];
  idx_t idx_shifted = idx - o2i[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(MAX(ic,idx_shifted),MIN(ic,idx_shifted),nc)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    } else cc[ic] = 0;
  }
}

static inline void // oc > ord
hpoly_der_gt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
//fprintf(stderr, "hpoly_der_gt: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[oc*ho + ord], *o2i = d->ord2idx;
  idx_t nc = o2i[oc+1] - o2i[oc];
  idx_t idx_shifted = idx - o2i[ord];
  for (idx_t ic = 0; ic < nc; ++ic) {
    idx_t ia = lc[hpoly_idx(idx_shifted,ic,nc)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
      *cnz = mad_bit_set(*cnz,oc);
    } else cc[ic] = 0;
  }
}

static inline void
hpoly_der(const T *a, idx_t idx, ord_t ord, T *c)
{
//fprintf(stderr, "hpoly_der: idx=%d, ord=%d\n", idx, ord);
  const D *d = c->d;
  const idx_t *o2i = d->ord2idx;
  const NUM *ca = a->coef;
        NUM *cc;

  c->hi = MIN3(c->mo, d->to, a->hi-ord);  // initial guess, readjust based on nz
  for (ord_t oc = 1; oc <= c->hi; ++oc) {
    if (mad_bit_tst(a->nz,oc+ord)) {
      cc = c->coef + o2i[oc];
           if (oc >  ord)  hpoly_der_gt(ca,cc,idx,oc,ord,&c->nz,d);
      else if (oc == ord)  hpoly_der_eq(ca,cc,idx,oc,ord,&c->nz,d);
      else   /*oc <  ord*/ hpoly_der_lt(ca,cc,idx,oc,ord,&c->nz,d);
    }
  }
}

// --- binary ops -------------------------------------------------------------o

// TPSA_LINOP(+, +, 0) => cc[i] = +ca[i] + cb[i], with i from (lo+0) to hi

#define TPSA_LINOP(OPA, OPB, ORD) \
do { \
  const idx_t *o2i = c->d->ord2idx; \
  idx_t start_a = o2i[ORD?MAX(a->lo,ORD):a->lo], end_a = o2i[MIN(a->hi,c_hi)+1]; \
  idx_t start_b = o2i[ORD?MAX(b->lo,ORD):b->lo], end_b = o2i[MIN(b->hi,c_hi)+1]; \
  idx_t i = start_a; \
  for (; i < MIN(end_a,start_b); ++i) c->coef[i] = OPA a->coef[i]; \
  for (; i <           start_b ; ++i) c->coef[i] = 0; \
  for (; i < MIN(end_a,end_b)  ; ++i) c->coef[i] = OPA a->coef[i] OPB b->coef[i]; \
  for (; i <     end_a         ; ++i) c->coef[i] = OPA a->coef[i]; \
  for (; i <           end_b   ; ++i) c->coef[i] =                OPB b->coef[i]; \
} while(0)

#define TPSA_CLEAR(ORD) \
do { \
  const idx_t *o2i = c->d->ord2idx; \
  for (idx_t i = o2i[ORD]; i < o2i[c_hi+1]; ++i) c->coef[i] = 0; \
} while(0)

void
FUN(scl) (const T *a, NUM v, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);

  const D *d = a->d;
  ensure(d == c->d, "incompatibles GTPSA (descriptors differ)");

  if (!v || a->hi == 0) {
    FUN(setvar)(c,v*a->coef[0],0,0); DBGFUN(<-); return;
  }

  FUN(copy0)(a,c);

  const idx_t *o2i = d->ord2idx;
  if (v == 1)
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
      c->coef[i] = a->coef[i];
  else if (v == -1)
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
      c->coef[i] = -a->coef[i];
  else
    for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
      c->coef[i] = v * a->coef[i];

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acc) (const T *a, NUM v, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);

  const D *d = c->d;
  ensure(a->d == d, "incompatibles GTPSA (descriptors differ)");

  if (!v) { DBGFUN(<-); return; }

  const idx_t *o2i = d->ord2idx;
  ord_t new_lo = MIN (a->lo,c->lo);
  ord_t new_hi = MIN3(a->hi,c->mo,d->to);

  c->nz = mad_bit_hcut(c->nz|a->nz, MAX(new_hi, c->hi));
  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  for (idx_t i = o2i[new_lo ]; i < o2i[c->lo   ]; ++i) c->coef[i]  = 0;
  for (idx_t i = o2i[c->hi+1]; i < o2i[new_hi+1]; ++i) c->coef[i]  = 0;
  for (idx_t i = o2i[a->lo  ]; i < o2i[new_hi+1]; ++i) c->coef[i] += v * a->coef[i];

  c->lo = new_lo;
  c->hi = MAX(new_hi, c->hi);
  if (c->lo) c->coef[0] = 0;
  if (TPSA_STRICT_NZ > 1) FUN(update0)(c, c->lo, c->hi);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(add) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(c);

  const D *d = a->d;
  ensure(d == b->d && d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, d->to);

  c->nz = mad_bit_hcut(a->nz|b->nz, c_hi);
  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }

  TPSA_LINOP(, +, 0);  // c->coef[i] = a->coef[i] + b->coef[i];

  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  if (c->lo) c->coef[0] = 0;
  if (TPSA_STRICT_NZ > 1) FUN(update0)(c, c->lo, c->hi);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sub) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(c);

  const D *d = a->d;
  ensure(d == b->d && d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t   hi = MAX (a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, d->to);

  c->nz = mad_bit_hcut(a->nz|b->nz, c_hi);
  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  const T* t = 0;
  if (a->lo > b->lo) SWAP(a,b,t);

  if (t) TPSA_LINOP(-, +, 0); // c->coef[i] = - a->coef[i] + b->coef[i];
  else   TPSA_LINOP( , -, 0); // c->coef[i] =   a->coef[i] - b->coef[i];

  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  if (c->lo) c->coef[0] = 0;
  if (TPSA_STRICT_NZ > 1) FUN(update0)(c, c->lo, c->hi);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(mul) (const T *a, const T *b, T *r)
{
  assert(a && b && r); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(r);

  const D *d = a->d;
  ensure(d == b->d && d == r->d, "incompatibles GTPSA (descriptors differ)");

  T *c = (a == r || b == r) ? GET_TMPX(r) : FUN(reset0)(r);

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
  if (!c->hi) goto ret;

  // order 1
  if (c->lo <= 1) {
    idx_t max_ord1 = d->ord2idx[2];
    if (mad_bit_tst(a->nz & b->nz,1) && a0 != 0 && b0 != 0) {
      for (idx_t i=1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i] + b0*a->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
    else if (mad_bit_tst(a->nz,1) && b0 != 0) {
      for (idx_t i=1; i < max_ord1; ++i) c->coef[i] = b0*a->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
    else if (mad_bit_tst(b->nz,1) && a0 != 0) {
      for (idx_t i=1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i];
      c->nz = mad_bit_set(c->nz,1);
    }
    else for (idx_t i=1; i < max_ord1; ++i) c->coef[i] = 0;
  }

  // order 2+
  if (c->hi > 1) {
    ord_t c_hi = c->hi, ab_hi = MAX(a->hi,b->hi);
    if (a0 != 0 && b0 != 0) TPSA_LINOP(b0*, +a0*, 2);
    else if       (b0 == 0) TPSA_LINOP( 0*, +a0*, 2);
    else                    TPSA_LINOP(b0*, + 0*, 2);
    TPSA_CLEAR(ab_hi+1);

    if (a0 != 0) { c->nz |= mad_bit_lcut(b->nz,2); }
    if (b0 != 0) { c->nz |= mad_bit_lcut(a->nz,2); }
                   c->nz  = mad_bit_hcut(c->nz,c_hi);

    if (mad_bit_tst(a->nz & b->nz,1)) {
      const idx_t *o2i = d->ord2idx, hod = d->mo/2;
      const idx_t *lc = d->L[hod+1];
      const idx_t *idx[2] = { d->L_idx[hod+1][0], d->L_idx[hod+1][2] };
      assert(lc);
      hpoly_diag_mul(a->coef+o2i[1], b->coef+o2i[1], c->coef, o2i[2]-o2i[1], lc, idx);
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

  if (TPSA_STRICT_NZ) FUN(update0)(c, c->lo, c->hi);

ret:
  assert(a != c && b != c);
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(div) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(c);
  ensure(a->d == b->d && a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  NUM b0 = b->coef[0];
  ensure(b0 != 0, "invalid domain");

  if (b->hi == 0) { FUN(scl)(a,1/b0,c); DBGFUN(<-); return; }

  T *t = (a == c || b == c) ? GET_TMPX(c) : FUN(reset0)(c);
  FUN(inv)(b,1,t);
  FUN(mul)(a,t,c);
  if (t != c) REL_TMPX(t);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(powi) (const T *a, int n, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);

  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  int inv = 0;

  if (n < 0) { n = -n; inv = 1; }

  T *t1 = GET_TMPX(c);

  switch (n) {
    case 0: FUN(setvar) (c, 1, 0, 0); break; // ok: no copy
    case 1: FUN(copy  ) (a, c);       break; // ok: 1 copy
    case 2: FUN(mul   ) (a,a, c);     break; // ok: 1 copy if a==c
    case 3: FUN(mul   ) (a,a, t1); FUN(mul)(t1,a,  c); break; // ok: 1 copy if a==c
    case 4: FUN(mul   ) (a,a, t1); FUN(mul)(t1,t1, c); break; // ok: no copy
    default: {
      T *t2 = GET_TMPX(c), *t;
      int ns = 0;
      FUN(copy  )(a, t1);
      FUN(setvar)(c, 1, 0, 0);
      for (;;) {
        if (n  & 1)   FUN(mul)(c ,t1, c ); // ok: 1 copy
        if (n /= 2) { FUN(mul)(t1,t1, t2); SWAP(t1,t2,t); ++ns; } // ok: no copy
        else break;
      }
      if (ns & 1) SWAP(t1,t2,t); // ensure even number of swaps
      REL_TMPX(t2);
    }
  }
  REL_TMPX(t1);

  if (inv) FUN(inv)(c,1,c);

  DBGTPSA(c); DBGFUN(<-);
}

log_t
FUN(equ) (const T *a, const T *b, num_t eps_)
{
  assert(a && b); DBGFUN(->); DBGTPSA(a); DBGTPSA(b);
  ensure(a->d == b->d, "incompatibles GTPSA (descriptors differ)");

  if (!(a->nz|b->nz)) return TRUE;

  if (eps_ < 0) eps_ = 1e-16;

  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }

  const idx_t *o2i = a->d->ord2idx;
  const idx_t start_a = o2i[a->lo], end_a = o2i[a->hi+1];
  const idx_t start_b = o2i[b->lo], end_b = o2i[b->hi+1];
  idx_t i = start_a;

  for (; i < MIN(end_a,start_b); ++i) if (fabs(a->coef[i]) > eps_)              { DBGFUN(<-); return FALSE; }
  i = start_b;
  for (; i < MIN(end_a,end_b)  ; ++i) if (fabs(a->coef[i] - b->coef[i]) > eps_) { DBGFUN(<-); return FALSE; }
  for (; i <     end_a         ; ++i) if (fabs(a->coef[i]) > eps_)              { DBGFUN(<-); return FALSE; }
  for (; i <           end_b   ; ++i) if (fabs(b->coef[i]) > eps_)              { DBGFUN(<-); return FALSE; }

  DBGFUN(<-); return TRUE;
}

// --- functions --------------------------------------------------------------o

#ifndef MAD_CTPSA_IMPL

void
FUN(abs) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  if (a->coef[0] < 0) FUN(scl) (a, -1, c);
  else if (a != c)    FUN(copy)(a, c);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(atan2) (const T *y, const T *x, T *r)
{
  assert(x && y && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(r);
  ensure(x->d == y->d && x->d == r->d, "incompatible GTPSA (descriptors differ)");
  NUM x0 = x->coef[0], y0 = y->coef[0];

  if (x0 != 0) {
    FUN(div)(y, x, r);
    FUN(atan)(r, r);
    if (x0 < 0) { // no copy ok
      mad_tpsa_scl(r, -1, r);
      if (y0 >= 0) FUN(set0)(r, 1,  M_PI_2);
      else         FUN(set0)(r, 1, -M_PI_2);
    }
  } else
    FUN(setvar)(r, atan2(y0,x0), 0,0); // Let C handle signs for ±pi/2

  DBGTPSA(r); DBGFUN(<-);
}

#else // MAD_CTPSA_IMPL

void
FUN(conj) (const T *a, T *c) // c = a.re - a.im I
{
  assert(a && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(c);
  ensure(a->d == c->d, "incompatibles GTPSA (descriptors differ)");

  c->lo = a->lo;
  c->hi = MIN3(a->hi, c->mo, c->d->to);
  c->nz = mad_bit_hcut(a->nz,c->hi);

  if (!c->nz) { FUN(reset0)(c); DBGFUN(<-); return; }

  const idx_t *o2i = c->d->ord2idx;
  for (idx_t i = o2i[c->lo]; i < o2i[c->hi+1]; ++i)
    c->coef[i] = conj(a->coef[i]);

  DBGTPSA(c); DBGFUN(<-);
}

#endif // MAD_CTPSA_IMPL

num_t
FUN(nrm) (const T *a)
{
  assert(a); DBGFUN(->); DBGTPSA(a);

  num_t nrm = 0;
  ord_t hi  = MIN(a->hi, a->d->to);
  if (mad_bit_hcut(a->nz,hi)) {
    const idx_t *o2i = a->d->ord2idx;
    for (idx_t i = o2i[a->lo]; i < o2i[hi+1]; ++i)
      nrm += fabs(a->coef[i]);
  }

  DBGFUN(<-); return nrm;
}

void
FUN(deriv) (const T *a, T *r, int iv)
{
  assert(a && r); DBGFUN(->); DBGTPSA(a); DBGTPSA(r);
  const D *d = a->d;
  const idx_t *o2i = d->ord2idx;

  ensure(d == r->d, "incompatibles GTPSA (descriptors differ)");
  ensure(iv >= o2i[1] && iv < o2i[2], "invalid domain");

  T *c = a == r ? GET_TMPX(r) : FUN(reset0)(r);

  if (!a->hi) goto ret; // empty

  FUN(setvar)(c,FUN(geti)(a,iv), 0, 0); // 0

  c->lo = a->lo ? a->lo-1 : 0;  // initial guess, readjusted after computation
  c->hi = MIN3(a->hi-1, c->mo, d->to);

  const NUM *ca = a->coef;

  ord_t der_ord = 1, oc = 1;
  if (mad_bit_tst(a->nz,oc+1))    // 1
      hpoly_der_eq(ca, c->coef+o2i[oc], iv, oc, der_ord, &c->nz, d);

  for (oc = 2; oc <= c->hi; ++oc) // 2..hi
    if (mad_bit_tst(a->nz,oc+1))
      hpoly_der_gt(ca, c->coef+o2i[oc], iv, oc, der_ord, &c->nz, d);

  if (c->nz) {
    c->lo = mad_bit_lowest (c->nz);
    c->hi = mad_bit_highest(c->nz);
    if (c->lo) c->coef[0] = 0;
    if (TPSA_STRICT_NZ) FUN(update0)(c, c->lo, c->hi);
  } else FUN(reset0)(c);

ret:
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(derivm) (const T *a, T *r, ssz_t n, const ord_t mono[n])
{
  assert(a && r); DBGFUN(->); DBGTPSA(a); DBGTPSA(r);
  const D *d = a->d;
  const idx_t *o2i = d->ord2idx;

  ensure(d == r->d, "incompatibles GTPSA (descriptors differ)");
  idx_t idx = mad_desc_idxm(d,n,mono);
  ensure(idx >= 0, "invalid monomial");

  // fallback on simple version
  if (idx < o2i[2]) { FUN(deriv)(a,r,idx); DBGFUN(<-); return; }

  T *c = a == r ? GET_TMPX(r) : FUN(reset0)(r);

  ord_t der_ord = mad_mono_ord(n,mono);
  ensure(der_ord > 0, "invalid derivative order");

  // ord 0 & setup
  FUN(setvar)(c,FUN(geti)(a,idx) * der_coef(idx,idx,der_ord,d), 0, 0);
  if (a->hi <= der_ord) goto ret;

  // ords 1..a->hi - 1
  hpoly_der(a, idx, der_ord, c);

  if (c->nz) {
    c->lo = mad_bit_lowest (c->nz);
    c->hi = mad_bit_highest(c->nz);
    if (c->lo) c->coef[0] = 0;
    if (TPSA_STRICT_NZ) FUN(update0)(c, c->lo, c->hi);
  } else FUN(reset0)(c);

ret:
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(poisson) (const T *a, const T *b, T *c, int nv)                 // C = [A,B]
{
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(c);
  const D *d = a->d;
  ensure(d == b->d && d == c->d, "incompatibles GTPSA (descriptors differ)");

  nv = nv>0 ? nv/2 : d->nv/2;

  T *is[3];
  for (int i=0; i < 3; ++i) is[i] = FUN(new)(a, d->to);

  FUN(reset0)(c);

  for (int i = 1; i <= nv; ++i) {
    FUN(deriv)(a, is[0], 2*i - 1); // res = res + da/dq_i * db/dp_i
    FUN(deriv)(b, is[1], 2*i    );
    FUN(mul)  (is[0],is[1],is[2]);
    FUN(add)  (c, is[2], c);

    FUN(deriv)(a, is[0], 2*i    ); // res = res - da/dp_i * db/dq_i
    FUN(deriv)(b, is[1], 2*i - 1);
    FUN(mul)  (is[0],is[1],is[2]);
    FUN(sub)  (c, is[2], c);
  }

  for (int i=0; i < 3; ++i) FUN(del)(is[i]);

  DBGTPSA(c); DBGFUN(<-);
}

// --- high level functions ---------------------------------------------------o

void
FUN(unit) (const T *x, T *r)
{
  assert(x && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(scl)(x, 1/fabs(x->coef[0]), r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(hypot) (const T *x, const T *y, T *r)
{
  assert(x && y && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(r);
  ensure(x->d == y->d && y->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(axypbvwpc)(1,x,x, 1,y,y, 0,r);
  FUN(sqrt)(r, r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(hypot3) (const T *x, const T *y, const T *z, T *r)
{
  assert(x && y && z && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(z); DBGTPSA(r);
  ensure(x->d == r->d && y->d == r->d && z->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(ax2pby2pcz2)(1,x, 1,y, 1,z, r);
  FUN(sqrt)(r, r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(axpb) (NUM a, const T *x, NUM b, T *r)
{
  assert(x && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(r);

  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");
  FUN(scl)(x,a,r);
  if (b) FUN(set0)(r,1,b);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(axpbypc) (NUM c1, const T *a, NUM c2, const T *b, NUM c3, T *c)
{            //    a           x       b           y      c      r
  assert(a && b && c); DBGFUN(->); DBGTPSA(a); DBGTPSA(b); DBGTPSA(c);
  ensure(a->d == b->d && b->d == c->d, "incompatibles GTPSA (descriptors differ)");

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->to);

  c->nz = mad_bit_hcut(a->nz|b->nz, c_hi);
  if (!c->nz) {
    FUN(setvar)(c, c1*a->coef[0]+c2*b->coef[0]+c3, 0,0); DBGFUN(<-); return;
  }

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

  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  if (c->lo) c->coef[0] = 0;
  if (c3) FUN(set0)(c,1,c3);
  if (TPSA_STRICT_NZ > 1) FUN(update0)(c, c->lo, c->hi);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(axypb) (NUM a, const T *x, const T *y, NUM b, T *r)
{
  assert(x && y && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(r);
  ensure(x->d == y->d && y->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(mul)(x,y,r);
  if (a != 1 || b != 0) FUN(axpb)(a,r,b,r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(axypbzpc) (NUM a, const T *x, const T *y, NUM b, const T *z, NUM c, T *r)
{
  assert(x && y && z && r);
  DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(z); DBGTPSA(r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = z == r ? GET_TMPX(r) : FUN(reset0)(r);
  FUN(mul)(x,y,t);
  FUN(axpbypc)(a,t,b,z,c,r);
  if (t != r) REL_TMPX(t);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(axypbvwpc) (NUM a, const T *x, const T *y,
                NUM b, const T *v, const T *w, NUM c, T *r)
{
  assert(x && y && v && w && r);
  DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(v); DBGTPSA(w); DBGTPSA(r);
  ensure(x->d == y->d && y->d == v->d && v->d == w->d && w->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(mul)(x,y,t);
  FUN(mul)(v,w,r);
  FUN(axpbypc)(a,t,b,r,c,r);
  REL_TMPX(t);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(ax2pby2pcz2) (NUM a, const T *x, NUM b, const T *y, NUM c, const T *z, T *r)
{
  assert(x && y && z && r);
  DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(z); DBGTPSA(r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d,
         "incompatibles GTPSA (descriptors differ)");

  T *t = z == r ? GET_TMPX(r) : FUN(reset0)(r);
  FUN(axypbvwpc)(a,x,x, b,y,y, 0, t);
  FUN(axypbzpc)(c,z,z, 1,t, 0, r);
  if (t != r) REL_TMPX(t);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(axpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");

  T *t = x == r ? GET_TMPX(r) : FUN(reset0)(r);
  FUN(axypb)(c,x,x,b,t);
  FUN(sqrt)(t,t);
  FUN(axpbypc)(a,x,1,t,0,r);
  if (t != r) REL_TMPX(t);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(logaxpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(r);
  ensure(x->d == r->d, "incompatibles GTPSA (descriptors differ)");

  FUN(axpsqrtbpcx2)(x, a, b, c, r);
  FUN(log)(r, r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(logxdy) (const T *x, const T *y, T *r)
{
  assert(x && y && r); DBGFUN(->); DBGTPSA(x); DBGTPSA(y); DBGTPSA(r);
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

  DBGTPSA(r); DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

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
