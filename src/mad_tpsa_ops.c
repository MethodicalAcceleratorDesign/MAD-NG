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

#include "mad_log.h"
#include "mad_cst.h"
#include "mad_num.h"
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
  FOR(ib,nb) if (cb[ib] || ca[ib])
    FOR(ia, idx[0][ib], idx[1][ib]) {
      idx_t ic = l[hpoly_idx(ib,ia,nb)];
      if (ic >= 0) cc[ic] += ca[ia]*cb[ib] + (ia != ib)*ca[ib]*cb[ia];
    }
}

static inline void
hpoly_sym_mul(const NUM *ca1, const NUM *cb1, const NUM *ca2, const NUM *cb2,
              NUM *cc, ssz_t na, ssz_t nb, const idx_t l[], const idx_t *idx[])
{
  // na > nb so longer loop is inside
  FOR(ib,nb) if (cb1[ib] || ca2[ib])
    FOR(ia, idx[0][ib], idx[1][ib]) {
      idx_t ic = l[hpoly_idx(ib,ia,na)];
      if (ic >= 0) cc[ic] += ca1[ia]*cb1[ib] + ca2[ib]*cb2[ia];
    }
}

static inline void
hpoly_asym_mul(const NUM *ca, const NUM *cb, NUM *cc, ssz_t na, ssz_t nb,
               const idx_t l[], const idx_t *idx[])
{
  // oa > ob so longer loop is inside
  FOR(ib,nb) if (cb[ib])
    FOR(ia, idx[0][ib], idx[1][ib]) {
      idx_t ic = l[hpoly_idx(ib,ia,na)];
      if (ic >= 0) cc[ic] += ca[ia]*cb[ib];
    }
}

static inline void
hpoly_mul(const T *a, const T *b, T *c, const ord_t *ocs, log_t in_parallel)
{
  const D *d = c->d;
  const idx_t *o2i = d->ord2idx;
  const NUM *ca = a->coef, *cb = b->coef;
  NUM   *cc = c->coef;
  idx_t hod = d->mo/2;
  bit_t nza = mad_bit_mask(~0ull, a->lo, a->hi);
  bit_t nzb = mad_bit_mask(~0ull, b->lo, b->hi);
  ord_t lo = MAX(c->lo,3);

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
        //printf("hpoly__sym_mul (%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_sym_mul(ca+o2i[oa],cb+o2i[ob],ca+o2i[ob],cb+o2i[oa],cc,na,nb,lc,idx);
      }
      else if (mad_bit_tst(nza,oa) && mad_bit_tst(nzb,ob)) {
        //printf("hpoly_asym_mul1(%d) %2d+%2d=%2d\n", ocs[0], oa,ob,oc);
        hpoly_asym_mul(ca+o2i[oa],cb+o2i[ob],cc,na,nb,lc,idx);
      }
      else if (mad_bit_tst(nza,ob) && mad_bit_tst(nzb,oa)) {
        //printf("hpoly_asym_mul2(%d) %2d+%2d=%2d\n", ocs[0], ob,oa,oc);
        hpoly_asym_mul(cb+o2i[oa],ca+o2i[ob],cc,na,nb,lc,idx);
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
      //printf("hpoly_diag_mul (%d) %2d+%2d=%2d\n", ocs[0], hoc,hoc,oc);
      hpoly_diag_mul(ca+o2i[hoc],cb+o2i[hoc],cc,nb,lc,idx);
    }
  }
}

#ifdef _OPENMP
static inline void
hpoly_mul_par(const T *a, const T *b, T *c) // parallel version
{
  const D *d = c->d;

  #pragma omp parallel for schedule(static,1)
  FOR(t,d->nth) {
    ord_t i = 0; while (d->ocs[1+t][i] > c->hi+1) ++i;
    hpoly_mul(a, b, c, &d->ocs[1+t][i], TRUE);
  }
}
#endif

static inline void
hpoly_mul_ser(const T *a, const T *b, T *c) // serial version
{
  hpoly_mul(a, b, c, &c->d->ocs[0][c->d->mo-c->hi], FALSE);
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
  if (mad_mono_lt(d->nv,srcm,derm)) return 0;

  num_t c = 1;
  FOR(v,d->nv)
    for (ord_t o = 0; o < derm[v]; ++o)
      c *= srcm[v] - o;
  return c;
}

static inline void // oc < ord
hpoly_der_lt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, const D *d)
{
//fprintf(stderr, "hpoly_der_lt: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *o2i = d->ord2idx;
  idx_t nc = o2i[oc+1] - o2i[oc], cols = o2i[ord+1] - o2i[ord];
  idx_t idx_shifted = idx - o2i[ord];
  FOR(ic,nc) {
    idx_t ia = lc[hpoly_idx(ic,idx_shifted,cols)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
    } else cc[ic] = 0;
  }
}

static inline void // oc == ord
hpoly_der_eq(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, const D *d)
{
//fprintf(stderr, "hpoly_der_eq: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *o2i = d->ord2idx;
  idx_t nc = o2i[ord+1] - o2i[ord];
  idx_t idx_shifted = idx - o2i[ord];
  FOR(ic,nc) {
    idx_t ia = lc[hpoly_idx(MAX(ic,idx_shifted),MIN(ic,idx_shifted),nc)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
    } else cc[ic] = 0;
  }
}

static inline void // oc > ord
hpoly_der_gt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, const D *d)
{
//fprintf(stderr, "hpoly_der_gt: idx=%d, oc=%d, ord=%d\n", idx, oc, ord);
  const idx_t ho = d->mo/2;
  const idx_t *lc = d->L[oc*ho + ord], *o2i = d->ord2idx;
  idx_t nc = o2i[oc+1] - o2i[oc];
  idx_t idx_shifted = idx - o2i[ord];
  FOR(ic,nc) {
    idx_t ia = lc[hpoly_idx(idx_shifted,ic,nc)];
    if (ia >= 0 && ca[ia]) {
      assert(o2i[oc+ord] <= ia && ia < o2i[oc+ord+1]);
      cc[ic] = ca[ia] * der_coef(ia,idx,ord,d);
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

  c->lo = a->lo > ord ? a->lo-ord : 1;
  c->hi = MIN(a->hi-ord, c->mo);
  if (c->lo > c->hi) { c->lo = 1, c->hi = 0; return; }

  for (ord_t oc = 1; oc <= c->hi; ++oc) {
    if (a->lo <= oc+ord && oc+ord <= a->hi) {
      cc = c->coef + o2i[oc];
           if (oc >  ord)  hpoly_der_gt(ca,cc,idx,oc,ord,d);
      else if (oc == ord)  hpoly_der_eq(ca,cc,idx,oc,ord,d);
      else   /*oc <  ord*/ hpoly_der_lt(ca,cc,idx,oc,ord,d);
    }
  }
}

// --- binary ops -------------------------------------------------------------o

static inline void
axpbypc (NUM c1, const T *a, NUM c2, const T *b, NUM c3, T *c)
{
  if (a == b) { FUN(scl)(a,c1+c2,c), c->val[0] += c3; return; }

  assert(a && b && c);

  log_t ac = FALSE, bc = FALSE;
  if (a == c) a = FUN(copy)(a,GETTMP(a)), ac = TRUE; else
  if (b == c) b = FUN(copy)(b,GETTMP(b)), bc = TRUE; else
  FUN(reset0)(c);

  c->val[0] = c1*a->val[0] + c2*b->val[0] + c3;

  idx_t ai=1, bi=1, ci=1;
  if (c1 && c2)
    while (ai < a->nc && bi < b->nc) {
      idx_t ii = a->idx[ai] - b->idx[bi];
      if (ii < 0) {
        if (a->idx[ai] >= c->mc) break;
        c->val[ci  ] = c1*a->val[ai  ];
        c->idx[ci++] =    a->idx[ai++];
      } else if (ii > 0) {
        if (b->idx[bi] >= c->mc) break;
        c->val[ci  ] = c2*b->val[bi  ];
        c->idx[ci++] =    b->idx[bi++];
      } else {
        if (a->idx[ai] >= c->mc) goto ret;
        NUM v = c1*a->val[ai] + c2*b->val[bi];
        if (v) c->val[ci] = v, c->idx[ci] = a->idx[ai], ci++;
        ai++, bi++;
      }
    }
  if (c1)
    while (ai < a->nc) {
      if (a->idx[ai] >= c->mc) break;
      c->val[ci  ] = c1*a->val[ai  ];
      c->idx[ci++] =    a->idx[ai++];
    }
  if (c2)
    while (bi < b->nc) {
      if (b->idx[bi] >= c->mc) break;
      c->val[ci  ] = c2*b->val[bi  ];
      c->idx[ci++] =    b->idx[bi++];
    }
ret:
  c->nc = ci;

  if (ac) RELTMP(a); else
  if (bc) RELTMP(b);
}

void
FUN(scl) (const T *a, NUM v, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  if (v == 0) { FUN(clear)(c);  DBGFUN(<-); return; }
  if (v == 1) { FUN(copy)(a,c); DBGFUN(<-); return; }

  if (FUN(isval)(a)) { FUN(setval)(c, v*a->val[0]); DBGFUN(<-); return; }

  c->nc = MIN(a->nc, c->mc);
  if (v == -1) { FOR(i,c->nc) c->val[i] =  -a->val[i]; }
  else         { FOR(i,c->nc) c->val[i] = v*a->val[i]; }

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(acc) (const T *a, NUM v, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  if (v == 0) { DBGFUN(<-); return; }

  axpbypc(v,a,1,c,0,c);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(add) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  axpbypc(1,a,1,b,0,c);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(sub) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  axpbypc(1,a,-1,b,0,c);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(dif) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  T *t = a == c ? GET_TMPX(c) : FUN(reset0)(c);

  axpbypc(1,a,-1,b,0,t);

  c->nc = MIN(t->nc, c->mc);
  idx_t nc = MIN(c->nc, a->nc);
  FOR(i,nc) c->val[i] = t->val[i] / MAX(fabs(a->val[i]),1);
  FOR(i,nc,c->nc) c->val[i] = t->val[i];

  if (t != c) REL_TMPX(t);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(axpbypc) (NUM c1, const T *a, NUM c2, const T *b, NUM c3, T *c)
{            //    a           x       b           y      c      r
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  axpbypc(c1,a,c2,b,c3,c);

  DBGTPSA(c); DBGFUN(<-);
}

void
FUN(mul) (const T *a, const T *b, T *r)
{
  assert(a && b && r); DBGFUN(->);
  ensure(IS_COMPAT(a,b,r), "incompatibles GTPSA (descriptors differ)");

  ord_t chi = MIN(a->hi+b->hi, r->mo);

  // order 0
  if (!chi) { FUN(setval)(r, a->coef[0]*b->coef[0]); DBGFUN(<-); return; }

  T *c = (a == r || b == r) ? GET_TMPX(r) : FUN(reset0)(r);

  if (a->lo > b->lo) { const T *t; SWAP(a,b,t); }

  // order 1+ and linear
  axpbypc(b->coef[0],a,a->coef[0],b,0,c), c->coef[0] *= 0.5;

  // order 2+
  if (chi > 1) {
    TPSA_SCAN(c,c->hi+1,chi) c->coef[i] = 0;
    c->lo = MIN(c->lo, a->lo+b->lo, c->mo);
    c->hi = chi;

    if (a->hi && b->hi && a->lo == 1 && b->lo == 1) {
      const D *d = a->d;
      const idx_t hod = d->mo/2;
      const idx_t *lc = d->L[hod+1];
      const idx_t *idx[2] = { d->L_idx[hod+1][0], d->L_idx[hod+1][2] };
      assert(lc);
      hpoly_diag_mul(a->coef+o2i[1], b->coef+o2i[1], c->coef, o2i[2]-o2i[1], lc, idx);
    }

    // order 3+
    if (chi > 2) {
#if !TPSA_STRICT
      FUN(nzero0)(a,a->lo,a->hi,1);
      FUN(nzero0)(b,b->lo,b->hi,1);
      if (a->lo > b->lo) { const T *t; SWAP(a,b,t); }
#endif

#ifdef _OPENMP
      if (o2i[c->hi+1] - o2i[c->lo] > 10000)
        hpoly_mul_par(a,b,c);
      else
#endif
        hpoly_mul_ser(a,b,c);

    }
#if TPSA_STRICT
  }
  FUN(update)(c);
#else
    FUN(update)(c);
  }
#endif

  assert(a != c && b != c);
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGFUN(<-);
}

void
FUN(div) (const T *a, const T *b, T *c)
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  NUM b0 = b->val[0];
  ensure(b0 != 0, "invalid domain");

  if (FUN(isval)(b)) {
#ifdef MAD_CTPSA_IMPL
    FUN(scl)(a, mad_cpx_inv(b0), c);
#else
    FUN(scl)(a, 1/b0, c);
#endif
    DBGFUN(<-); return;
  }

  T *t = (a == c || b == c) ? GET_TMPX(c) : FUN(reset0)(c);
  FUN(inv)(b,1,t);
  FUN(mul)(a,t,c);
  if (t != c) REL_TMPX(t);
  DBGFUN(<-);
}

void
FUN(powi) (const T *a, int n, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  int inv = 0;

  if (n < 0) { n = -n; inv = 1; }

  T *t1 = n > 2 ? GET_TMPX(c) : NULL;

  switch (n) {
    case 0: FUN(setval)(c,  1);  break; // ok: no copy
    case 1: FUN(copy  )(a,  c);  break; // ok: 1 copy if a!=c
    case 2: FUN(mul   )(a,a,c);  break; // ok: 1 copy if a==c
    case 3: FUN(mul   )(a,a,t1); FUN(mul)(t1,a, c); break; // ok: 1 copy if a==c
    case 4: FUN(mul   )(a,a,t1); FUN(mul)(t1,t1,c); break; // ok: no copy
    default: {
      T *t2 = GET_TMPX(c), *t;
      int ns = 0;
      FUN(copy)(a, t1);
      FUN(setval)(c, 1);
      for (;;) {
        if (n  & 1)   FUN(mul)(c ,t1,c ); // ok: 1 copy
        if (n /= 2) { FUN(mul)(t1,t1,t2); SWAP(t1,t2,t); ++ns; } // ok: no copy
        else break;
      }
      if (ns & 1) SWAP(t1,t2,t); // ensure even number of swaps
      REL_TMPX(t2);
    }
  }
  if (t1) REL_TMPX(t1);

  if (inv) FUN(inv)(c,1,c);
  DBGFUN(<-);
}

static inline log_t
neq (NUM a, NUM b, num_t tol) {
  NUM d = a - b;
#ifndef MAD_CTPSA_IMPL
  return fabs(d) > tol; }
#else
  return fabs(creal(d)) > tol || fabs(cimag(d)) > tol; }
#endif

log_t
FUN(equ) (const T *a, const T *b, num_t tol)
{
  assert(a && b); DBGFUN(->);
  ensure(IS_COMPAT(a,b), "incompatibles GTPSA (descriptors differ)");
  ensure(FALSE, "NYI");

  const D *d = a->d;
  T c_ = {.d=d, .mo=d->mo, .ao=d->mo}, *c = &c_; // fake TPSA

  if (a->lo > b->lo) { const T *t; SWAP(a,b,t); }

  ord_t alo = a->lo, ahi = a->hi;
  ord_t blo = b->lo, bhi = b->hi;
  FUN(copy00)(a,b,c);

  log_t equ = !neq(a->coef[0], b->coef[0], tol);

  if (!c->hi || !equ) { DBGFUN(<-); return equ; }

  if (ahi > c->hi) ahi = c->hi;
  if (bhi > c->hi) bhi = c->hi;

  const idx_t *o2i = a->d->ord2idx;
  idx_t i = o2i[alo];
  for(; i < o2i[MIN(ahi+1,blo)]; i++) if (neq(a->coef[i],0         ,tol)) goto ret;
  for(; i < o2i[MIN(ahi,bhi)+1]; i++) if (neq(a->coef[i],b->coef[i],tol)) goto ret;
  for(; i < o2i[ahi+1]         ; i++) if (neq(a->coef[i],0         ,tol)) goto ret;
  for(; i < o2i[bhi+1]         ; i++) if (neq(0         ,b->coef[i],tol)) goto ret;

  DBGFUN(<-); return TRUE;
ret:
  DBGFUN(<-); return FALSE;
}

// --- functions --------------------------------------------------------------o

#ifndef MAD_CTPSA_IMPL

void
FUN(abs) (const T *a, T *c)
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  if (a->val[0] < 0) FUN(scl) (a, -1, c);
  else if (a != c)   FUN(copy)(a, c);

  DBGFUN(<-);
}

void
FUN(atan2) (const T *y, const T *x, T *r)
{
  assert(x && y && r); DBGFUN(->);
  ensure(IS_COMPAT(y,x,r), "incompatibles GTPSA (descriptors differ)");

  NUM x0 = x->val[0], y0 = y->val[0];
  NUM a0 = atan2(y0, x0);

  if (fabs(a0) < M_PI_2-0.1 || fabs(a0) > M_PI_2+0.1) {
    FUN(div)(y, x, r);
    FUN(atan)(r, r);
  } else {
    FUN(axypbvwpc)(1,x,x, 1,y,y, 0, r);
    FUN(invsqrt)(r, 1, r);
    FUN(mul)(x, r, r);
    FUN(acos)(r, r);
    if (y0 < 0) FUN(scl)(r, -1, r);
  }
  FUN(seti)(r, 0, 0, a0);

  DBGFUN(<-);
}

#else // MAD_CTPSA_IMPL

void
FUN(conj) (const T *a, T *c) // c = a.re - a.im I
{
  assert(a && c); DBGFUN(->);
  ensure(IS_COMPAT(a,c), "incompatibles GTPSA (descriptors differ)");

  c->nc = MIN(a->nc, c->mc);
  FOR(i,c->nc) c->val[i] = conj(a->val[i]);

  DBGTPSA(c); DBGFUN(<-);
}

#endif // MAD_CTPSA_IMPL

num_t
FUN(nrm) (const T *a)
{
  assert(a); DBGFUN(->);
  num_t nrm = 0;
  FOR(i,a->nc) nrm += fabs(a->val[i]);
  DBGFUN(<-); return nrm;
}

void
FUN(integ) (const T *a, T *r, int iv)
{
  assert(a && r); DBGFUN(->);
  const D *d = a->d;
  ensure(IS_COMPAT(a,r), "incompatibles GTPSA (descriptors differ)");
  ensure(0 < iv && iv <= d->nv,
         "index 1<= %d <=%d is not a GTPSA variable", iv, d->nv);

  T *c = a == r ? GET_TMPX(r) : FUN(reset0)(r);

  T *t = GET_TMPX(r);
  FUN(setvar)(t, 0, iv, 0);
  FUN(mul)(a, t, c);    // integrate
  REL_TMPX(t);

  ord_t **mono = d->To;
  TPSA_SCAN(c,MAX(c->lo,2),c->hi)
    if (c->coef[i] && mono[i][iv-1] > 1) c->coef[i] /= mono[i][iv-1];

  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGFUN(<-);
}

void
FUN(deriv) (const T *a, T *r, int iv)
{
  assert(a && r); DBGFUN(->);
  const D *d = a->d;
  ensure(IS_COMPAT(a,r), "incompatibles GTPSA (descriptors differ)");
  ensure(0 < iv && iv <= d->nv,
         "index 1<= %d <=%d is not a GTPSA variable", iv, d->nv);

  if (FUN(isval)(a)) { FUN(clear)(r); DBGFUN(<-); return; } // empty

  T *c = a == r ? GET_TMPX(r) : FUN(reset0)(r);

  FUN(setval)(c, FUN(geti)(a,iv));

  c->lo = a->lo > 1 ? a->lo-1 : 1;
  c->hi = MIN(a->hi-1, c->mo);
  if (c->lo > c->hi) { c->lo = 1, c->hi = 0; goto ret; }

  const idx_t *o2i = d->ord2idx;
  ord_t der_ord = 1, oc = 1;
  if (a->lo <= oc+1 && oc+1 <= a->hi) // 1
      hpoly_der_eq(a->coef, c->coef+o2i[oc], iv, oc, der_ord, d);

  for (oc = 2; oc <= c->hi; ++oc) // 2..hi
    if (a->lo <= oc+1 && oc+1 <= a->hi)
      hpoly_der_gt(a->coef, c->coef+o2i[oc], iv, oc, der_ord, d);

  FUN(update)(c);

ret:
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGFUN(<-);
}

void
FUN(derivm) (const T *a, T *r, ssz_t n, const ord_t mono[n])
{
  assert(a && r); DBGFUN(->);
  ensure(IS_COMPAT(a,r), "incompatibles GTPSA (descriptors differ)");

  const D *d = a->d;
  idx_t idx = mad_desc_idxm(d,n,mono);
  ensure(idx >= 0, "invalid monomial");

  // fallback on simple version
  if (idx < d->ord2idx[2]) { FUN(deriv)(a,r,idx); DBGFUN(<-); return; }

  T *c = a == r ? GET_TMPX(r) : FUN(reset0)(r);

  ord_t der_ord = mad_mono_ord(n,mono);
  ensure(der_ord > 0, "invalid derivative order");

  // ord 0 & setup
  FUN(setval)(c, FUN(geti)(a,idx) * der_coef(idx,idx,der_ord,d));
  if (a->hi <= der_ord) goto ret;

  // ords 1..a->hi - 1
  hpoly_der(a, idx, der_ord, c);

  FUN(update)(c);

ret:
  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGFUN(<-);
}

void
FUN(poisbra) (const T *a, const T *b, T *r, int nv)                 // C = [A,B]
{
  assert(a && b && r); DBGFUN(->); DBGTPSA(a); DBGTPSA(b);
  ensure(IS_COMPAT(a,b,r), "incompatibles GTPSA (descriptors differ)");

  nv = nv>0 ? nv/2 : a->d->nv/2;

  T *c = a == r || b == r ? GET_TMPX(r) : FUN(reset0)(r);

  T *is[3];
  FOR(i,3) is[i] = FUN(new)(a, mad_tpsa_same);
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
  FOR(i,3) FUN(del)(is[i]);

  if (c != r) { FUN(copy)(c,r); REL_TMPX(c); }
  DBGFUN(<-);
}

// --- high level functions ---------------------------------------------------o

void
FUN(unit) (const T *x, T *r)
{
  assert(x && r); DBGFUN(->);
  ensure(IS_COMPAT(x,r), "incompatibles GTPSA (descriptors differ)");

  FUN(scl)(x, 1/fabs(x->coef[0]), r);

  DBGFUN(<-);
}

void
FUN(hypot) (const T *x, const T *y, T *r)
{
  assert(x && y && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,r), "incompatibles GTPSA (descriptors differ)");

  FUN(axypbvwpc)(1,x,x, 1,y,y, 0,r);
  FUN(sqrt)(r, r);

  DBGFUN(<-);
}

void
FUN(hypot3) (const T *x, const T *y, const T *z, T *r)
{
  assert(x && y && z && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,z,r), "incompatibles GTPSA (descriptors differ)");

  FUN(ax2pby2pcz2)(1,x, 1,y, 1,z, r);
  FUN(sqrt)(r, r);

  DBGFUN(<-);
}

void
FUN(axpb) (NUM a, const T *x, NUM b, T *r)
{
  assert(x && r); DBGFUN(->);
  ensure(IS_COMPAT(x,r), "incompatibles GTPSA (descriptors differ)");

  FUN(scl)(x,a,r);
  if (b) FUN(seti)(r,0,1,b);

  DBGFUN(<-);
}

void
FUN(axypb) (NUM a, const T *x, const T *y, NUM b, T *r)
{
  assert(x && y && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,r), "incompatibles GTPSA (descriptors differ)");

  FUN(mul)(x,y,r);
  if (a != 1 || b != 0) FUN(axpb)(a,r,b,r);

  DBGFUN(<-);
}

void
FUN(axypbzpc) (NUM a, const T *x, const T *y, NUM b, const T *z, NUM c, T *r)
{
  assert(x && y && z && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,z,r), "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(mul)(x,y,t);
  FUN(axpbypc)(a,t,b,z,c,r);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(axypbvwpc) (NUM a, const T *x, const T *y,
                NUM b, const T *v, const T *w, NUM c, T *r)
{
  assert(x && y && v && w && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,v,w,r), "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(mul)(x,y,t);
  FUN(mul)(v,w,r);
  FUN(axpbypc)(a,t,b,r,c,r);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(ax2pby2pcz2) (NUM a, const T *x, NUM b, const T *y, NUM c, const T *z, T *r)
{
  assert(x && y && z && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,z,r), "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(axypbvwpc)(a,x,x,b,y,y,0,t);
  FUN(axypbzpc)(c,z,z,1,t,0,r);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(axpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r); DBGFUN(->);
  ensure(IS_COMPAT(x,r), "incompatibles GTPSA (descriptors differ)");

  T *t = GET_TMPX(r);
  FUN(axypb)(c,x,x,b,t);
  FUN(sqrt)(t,t);
  FUN(axpbypc)(a,x,1,t,0,r);
  REL_TMPX(t); DBGFUN(<-);
}

void
FUN(logaxpsqrtbpcx2) (const T *x, NUM a, NUM b, NUM c, T *r)
{
  assert(x && r); DBGFUN(->);
  ensure(IS_COMPAT(x,r), "incompatibles GTPSA (descriptors differ)");

  FUN(axpsqrtbpcx2)(x,a,b,c,r);
  FUN(log)(r,r);
  DBGFUN(<-);
}

void
FUN(logxdy) (const T *x, const T *y, T *r)
{
  assert(x && y && r); DBGFUN(->);
  ensure(IS_COMPAT(x,y,r), "incompatibles GTPSA (descriptors differ)");

  NUM x0 = FUN(geti)(x,0), y0 = FUN(geti)(y,0);
  if (fabs(x0) > fabs(y0)) { // TODO: improve stability around 1i
    FUN(div)(x, y, r);
    FUN(log)(r, r);
  } else {
    FUN(div)(y, x, r);
    FUN(log)(r, r);
    FUN(scl)(r, -1, r);
  }
  DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(acc_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(acc)(a, CPX(v), c); }

void FUN(scl_r) (const T *a, num_t v_re, num_t v_im, T *c)
{ FUN(scl)(a, CPX(v), c); }

void FUN(axpb_r) (num_t a_re, num_t a_im, const T *x,
                  num_t b_re, num_t b_im, T *r)
{ FUN(axpb)(CPX(a), x, CPX(b), r); }

void FUN(axpbypc_r) (num_t a_re, num_t a_im, const T *x,
                     num_t b_re, num_t b_im, const T *y,
                     num_t c_re, num_t c_im, T *r)
{ FUN(axpbypc)(CPX(a), x, CPX(b), y, CPX(c), r); }

void FUN(axypb_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                   num_t b_re, num_t b_im, T *r)
{ FUN(axypb)(CPX(a), x, y, CPX(b), r); }

void FUN(axypbzpc_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                      num_t b_re, num_t b_im, const T *z,
                      num_t c_re, num_t c_im, T *r)
{ FUN(axypbzpc)(CPX(a), x, y, CPX(b), z, CPX(c), r); }

void FUN(axypbvwpc_r) (num_t a_re, num_t a_im, const T *x, const T *y,
                       num_t b_re, num_t b_im, const T *v, const T *w,
                       num_t c_re, num_t c_im, T *r)
{ FUN(axypbvwpc)(CPX(a), x, y, CPX(b), v, w, CPX(c), r); }

void FUN(ax2pby2pcz2_r)(num_t a_re, num_t a_im, const T *x,
                        num_t b_re, num_t b_im, const T *y,
                        num_t c_re, num_t c_im, const T *z, T *r)
{ FUN(ax2pby2pcz2)(CPX(a), x, CPX(b), y, CPX(c), z, r); }

void FUN(axpsqrtbpcx2_r)(const T *x, num_t a_re, num_t a_im,
                                     num_t b_re, num_t b_im,
                                     num_t c_re, num_t c_im, T *r)
{ FUN(axpsqrtbpcx2)(x, CPX(a), CPX(b), CPX(c), r); }

void FUN(logaxpsqrtbpcx2_r) (const T *x, num_t a_re, num_t a_im,
                                         num_t b_re, num_t b_im,
                                         num_t c_re, num_t c_im, T *r)
{ FUN(logaxpsqrtbpcx2)(x, CPX(a), CPX(b), CPX(c), r); }

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
