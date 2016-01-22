/*
 o----------------------------------------------------------------------------o
 |
 | TPSA operators module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
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

// --- LOCAL FUNCTIONS --------------------------------------------------------

static inline void
hpoly_triang_mul(const NUM *ca, const NUM *cb, NUM *cc, int nb,
                 const idx_t l[], const int *idx[])
{
  // asymm: c[2 2] = a[2 0]*b[0 2] + a[0 2]*b[2 0]
  for (idx_t ib = 0; ib < nb; ib++)
    if (cb[ib] || ca[ib])
      for (idx_t ia = idx[0][ib]; ia < idx[1][ib]; ia++) {
        int ic = l[hpoly_idx(ib,ia,nb)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca[ia]*cb[ib] + (ia == ib ? 0 : ca[ib]*cb[ia]);
      }
}

static inline void
hpoly_sym_mul(const NUM *ca1, const NUM *cb1, const NUM *ca2, const NUM *cb2,
              NUM *cc, int na, int nb, const idx_t l[], const int *idx[])
{
  // na > nb so longer loop is inside
  for (idx_t ib=0; ib < nb; ib++)
    if (cb1[ib] || ca2[ib])
      for (idx_t ia=idx[0][ib]; ia < idx[1][ib]; ia++) {
        int ic = l[hpoly_idx(ib, ia, na)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca1[ia]*cb1[ib] + ca2[ib]*cb2[ia];
      }
}

static inline void
hpoly_asym_mul(const NUM *ca, const NUM *cb, NUM *cc, int na, int nb,
               const idx_t l[], const int *idx[])
{
  // oa > ob so longer loop is inside
  for (idx_t ib=0; ib < nb; ib++)
    if (cb[ib])
      for (idx_t ia=idx[0][ib]; ia < idx[1][ib]; ia++) {
        int ic = l[hpoly_idx(ib,ia,na)];
        if (ic >= 0)
          cc[ic] = cc[ic] + ca[ia]*cb[ib];
      }
}

static inline void
hpoly_mul(const T *a, const T *b, T *c, const ord_t *ocs, bit_t *cnz, int in_parallel)
{
  D *d = c->d;
  int *pi = d->ord2idx, hod = d->mo/2;
  const NUM *ca = a->coef,  *cb = b->coef;
        NUM *cc = c->coef;
        bit_t nza = a->nz  ,  nzb = b->nz;

  for (int i = 0; ocs[i]; ++i) {
    if (ocs[i] < c->lo || ocs[i] > c->hi + 1 || (ocs[i] == c->hi + 1 && !in_parallel))
      continue;

    ord_t oc = ocs[i];
    int idx0 = 0, idx1 = 2;
    if (in_parallel && ocs[i] >= c->hi) {
      oc = c->hi;
      if (ocs[i] == c->hi) idx1 = 1;
      else                 idx0 = 1;
    }

    for (int j=1; j <= (oc-1)/2; ++j) {
      int oa = oc-j, ob = j;            // oa > ob >= 1
      int na = pi[oa+1] - pi[oa];
      int nb = pi[ob+1] - pi[ob];
      const idx_t *lc  = d->L[oa*hod + ob];
      assert(lc);
      const int *idx[2] = { d->L_idx[oa*hod + ob][idx0],
                            d->L_idx[oa*hod + ob][idx1]};
      assert(idx[0] && idx[1]);

      if (mad_bit_get(nza,oa) && mad_bit_get(nzb,ob) &&
          mad_bit_get(nza,ob) && mad_bit_get(nzb,oa)) {
        hpoly_sym_mul(ca+pi[oa],cb+pi[ob], ca+pi[ob],cb+pi[oa], cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_get(nza,oa) && mad_bit_get(nzb,ob)) {
        hpoly_asym_mul(ca+pi[oa],cb+pi[ob],cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
      else if (mad_bit_get(nza,ob) && mad_bit_get(nzb,oa)) {
        hpoly_asym_mul(cb+pi[oa],ca+pi[ob],cc, na,nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
    }

    if (! (oc & 1)) {  // even oc, triang matrix
      int hoc = oc/2, nb = pi[hoc+1]-pi[hoc];
      const idx_t *lc = d->L[hoc*hod + hoc];
      const int *idx[2] = { d->L_idx[hoc*hod + hoc][idx0],
                            d->L_idx[hoc*hod + hoc][idx1] };
      assert(lc);
      if (mad_bit_get(nza,hoc) && mad_bit_get(nzb,hoc) ) {
        hpoly_triang_mul(ca+pi[hoc],cb+pi[hoc],cc, nb, lc, idx);
        *cnz = mad_bit_set(*cnz,oc);
      }
    }
  }
}

#ifdef _OPENMP
static inline void
hpoly_mul_par(const T *a, const T *b, T *c)
{
  int nb_threads = omp_get_num_procs();
  bit_t c_nzs[nb_threads];
  for (int t = 0; t < nb_threads; ++t)
    c_nzs[t] = c->nz;

  #pragma omp parallel for
  for (int t = 0; t < nb_threads; ++t)
    hpoly_mul(a,b,c,c->d->ocs[t],&c_nzs[t],1);

  for (int t = 0; t < nb_threads; ++t)
    c->nz |= c_nzs[t];
}
#endif

static inline void
hpoly_mul_ser(const T *a, const T *b, T *c)
{
  hpoly_mul(a,b,c,c->d->ocs[0],&c->nz,0);
}

static inline int
der_coef(idx_t ia, idx_t di, ord_t der_ord, const D* d)
{
  if (der_ord == 1) // der against var
    return d->To[ia][di-1];

  const ord_t *msrc = d->To[ia], *mder = d->To[di];
  if (!mad_mono_leq(d->nv,mder,msrc))
    return 0;

  int c = 1;
  for (int v = 0; v < d->nv; ++v)
    for (int o = 0; o < mder[v]; ++o)
        c *= msrc[v] - o;
  return c;
}

static inline void
hpoly_der_lt(const NUM ca[], NUM cc[], idx_t idx, ord_t oc, ord_t ord, bit_t *cnz, const D *d)
{
  const ord_t ho = d->mo/2;
  const idx_t *lc = d->L[ord*ho + oc], *pi = d->ord2idx;
  int nc = pi[oc+1] - pi[oc], cols = pi[ord+1] - pi[ord];
  // idx = idx - pi[ord];
  for (int ic = 0; ic < nc; ++ic) {
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
  int nc = pi[ord+1] - pi[ord];
  idx_t idx_shifted = idx - pi[ord];
  for (int ic = 0; ic < nc; ++ic) {
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
  int nc = pi[oc+1] - pi[oc];
  idx_t idx_shifted = idx - pi[ord];
  for (int ic = 0; ic < nc; ++ic) {
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
  D *d = c->d;
  idx_t *pi = d->ord2idx;
  const NUM *ca = a->coef;
        NUM *cc;

  c->hi = MIN3(c->mo, d->trunc, a->hi-ord);  // initial guess, readjust based on nz
  for (int oc = 1; oc <= c->hi; ++oc)
    if (mad_bit_get(a->nz,oc+ord)) {
      cc = c->coef + pi[oc];
      if (oc < ord)
        hpoly_der_lt(ca,cc,idx,oc,ord,&c->nz,d);
      else if (oc == ord)
        hpoly_der_eq(ca,cc,idx,oc,ord,&c->nz,d);
      else
        hpoly_der_gt(ca,cc,idx,oc,ord,&c->nz,d);
    }
  int n = mad_bit_lowest(c->nz);
  c->lo = MIN(n,c->mo);
  c->hi = mad_bit_highest(c->nz);
}

// --- PUBLIC FUNCTIONS -------------------------------------------------------

// --- --- UNARY ---------------------------------------------------------------

void
FUN(abs) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);

  c->hi = MIN3(a->hi, c->mo, c->d->trunc);
  c->lo = a->lo;
  c->nz = mad_bit_trunc(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i) {
    c->coef[i] = fabs(a->coef[i]);
  }
}

#ifdef MAD_CTPSA_IMPL

void
FUN(arg) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);

  c->hi = MIN3(a->hi, c->mo, c->d->trunc);
  c->lo = a->lo;
  c->nz = mad_bit_trunc(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = carg(a->coef[i]);
}

void
FUN(conj) (const T *a, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);

  c->hi = MIN3(a->hi, c->mo, c->d->trunc);
  c->lo = a->lo;
  c->nz = mad_bit_trunc(a->nz,c->hi);

  idx_t *pi = c->d->ord2idx;
  for (int i = pi[c->lo]; i < pi[c->hi+1]; ++i)
    c->coef[i] = conj(a->coef[i]);
}

#endif

NUM
FUN(nrm1) (const T *a, const T *b_)
{
  assert(a);
  NUM norm = 0.0;
  int *pi = a->d->ord2idx;
  if (b_) {
    ensure(a->d == b_->d);
    if (a->lo > b_->lo) { const T *t; SWAP(a,b_,t); }

    idx_t start_a = pi[a ->lo], end_a = pi[MIN(a ->hi,a ->d->trunc)+1],
          start_b = pi[b_->lo], end_b = pi[MIN(b_->hi,b_->d->trunc)+1];
    idx_t i;
    for (i = start_a; i < MIN(end_a,start_b); ++i) norm += fabs(a->coef[i]);
    for (i = start_b; i < MIN(end_a,end_b)  ; ++i) norm += fabs(a->coef[i] - b_->coef[i]);
    for (           ; i <     end_a         ; ++i) norm += fabs(a->coef[i]);
    for (           ; i <     end_b         ; ++i) norm += fabs(b_->coef[i]);
  }
  else {
    ord_t hi = MIN(a->hi, a->d->trunc);
    for (int o = a->lo; o <= hi; ++o)
      if (mad_bit_get(a->nz,o)) {
        for (int i = pi[o]; i < pi[o+1]; ++i)
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
  int *pi = a->d->ord2idx;
  if (b_) {
    ensure(a->d == b_->d);
    if (a->lo > b_->lo) { const T* t; SWAP(a,b_,t); }

    idx_t start_a = pi[a ->lo], end_a = pi[MIN(a ->hi,a ->d->trunc)+1],
          start_b = pi[b_->lo], end_b = pi[MIN(b_->hi,b_->d->trunc)+1];
    idx_t i;
    for (i = start_a; i < MIN(end_a,start_b); ++i) norm +=  a->coef[i]              *  a->coef[i];
    for (i = start_b; i < MIN(end_a,end_b)  ; ++i) norm += (a->coef[i]-b_->coef[i]) * (a->coef[i]-b_->coef[i]);
    for (           ; i <     end_a         ; ++i) norm +=  a->coef[i]              *  a->coef[i];
    for (           ; i <     end_b         ; ++i) norm +=  b_->coef[i]             *  b_->coef[i];
  }
  else {
    ord_t hi = MIN(a->hi, a->d->trunc);
    for (int o = a->lo; o <= hi; ++o)
      if (mad_bit_get(a->nz,o)) {
        for (int i = pi[o]; i < pi[o+1]; ++i)
          norm += a->coef[i] * a->coef[i];
      }
  }
  return sqrt(norm);
}

void
FUN(der) (const T *a, T *c, int var)
{
  assert(a && c && a != c);
  ensure(var >= a->d->ord2idx[1] && var < a->d->ord2idx[2]);
  // TODO: ensure map_order[var] > 0

  if (a->hi == 0) { FUN(clear)(c); return; }
  FUN(scalar)(c,FUN(geti)(a,var));  // TODO: what if alpha[var] == 0 ?

  D *d = c->d;
  c->hi = MIN3(c->mo, d->trunc, a->hi-1);
  c->lo = a->lo ? a->lo-1 : 0;  // initial guess, readjusted after computation

  idx_t *pi = d->ord2idx;
  const NUM *ca = a->coef;

  ord_t der_ord = 1, oc = 1;
  if (mad_bit_get(a->nz,oc+1))
      hpoly_der_eq(ca,c->coef+pi[oc],var,oc,der_ord,&c->nz,d);
  for (oc = 2; oc <= c->hi; ++oc)
    if (mad_bit_get(a->nz,oc+1))
      hpoly_der_gt(ca,c->coef+pi[oc],var,oc,der_ord,&c->nz,d);
  int n = mad_bit_lowest(c->nz);
  c->lo = MIN(n,c->mo);
  c->hi = mad_bit_highest(c->nz);
}

void
FUN(mder) (const T *a, T *c, int n, const ord_t mono[n])
{
  assert(a && c && a != c);
  assert(a->d == c->d);
  ensure(mad_desc_mono_isvalid(a->d,n,mono));

  ord_t der_ord = mad_mono_ord(n,mono);
  ensure(der_ord > 0);
  idx_t idx = mad_desc_get_idx(a->d,n,mono);
  if (idx < a->d->ord2idx[2]) {  // fallback on simple version
    FUN(der)(a,c,idx);
    return;
  }

  // ord 0 & setup
  FUN(scalar)(c,FUN(geti)(a,idx) * der_coef(idx,idx,der_ord,a->d));
  if (a->hi <= der_ord)
    return;

  // ords 1..a->hi - 1
  hpoly_der(a,idx,der_ord,c);
}

void
FUN(scl) (const T *a, NUM v, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);

  if (a->hi == 0) { FUN(scalar)(c, v*a->coef[0]); return; }

  D *d = a->d;
  c->lo = a->lo;
  c->hi = MIN3(a->hi, c->mo, d->trunc);
  c->nz = mad_bit_trunc(a->nz,c->hi);
  for (int i = d->ord2idx[c->lo]; i < d->ord2idx[c->hi+1]; ++i)
    c->coef[i] = v * a->coef[i];
}

// --- --- BINARY --------------------------------------------------------------

// TPSA_LINOP_ORD(+, +, 0) => cc[i] = +ca[i] + cb[i], with i from (lo+0) to hi
// TPSA_LINOP assumes ORD=0 and avoids GCC warning -Wtype-limits

#define TPSA_LINOP(OPA, OPB) \
do { \
    idx_t *pi = c->d->ord2idx; \
    idx_t start_a = pi[a->lo], end_a = pi[MIN(a->hi,c_hi)+1]; \
    idx_t start_b = pi[b->lo], end_b = pi[MIN(b->hi,c_hi)+1]; \
    int i = start_a; \
    for (; i < MIN(end_a,start_b); ++i) c->coef[i] = OPA a->coef[i]; \
    for (; i <           start_b ; ++i) c->coef[i] = 0; \
    for (; i < MIN(end_a,end_b)  ; ++i) c->coef[i] = OPA a->coef[i] OPB b->coef[i]; \
    for (; i <     end_a         ; ++i) c->coef[i] = OPA a->coef[i]; \
    for (; i <           end_b   ; ++i) c->coef[i] =                OPB b->coef[i]; \
} while(0) \

#define TPSA_LINOP_ORD(OPA, OPB, ORD) \
do { \
    idx_t *pi = c->d->ord2idx; \
    idx_t start_a = pi[MAX(a->lo,ORD)], end_a = pi[MIN(a->hi,c_hi)+1]; \
    idx_t start_b = pi[MAX(b->lo,ORD)], end_b = pi[MIN(b->hi,c_hi)+1]; \
    int i = start_a; \
    for (; i < MIN(end_a,start_b); ++i) c->coef[i] = OPA a->coef[i]; \
    for (; i <           start_b ; ++i) c->coef[i] = 0; \
    for (; i < MIN(end_a,end_b)  ; ++i) c->coef[i] = OPA a->coef[i] OPB b->coef[i]; \
    for (; i <     end_a         ; ++i) c->coef[i] = OPA a->coef[i]; \
    for (; i <           end_b   ; ++i) c->coef[i] =                OPB b->coef[i]; \
} while(0) \

void
FUN(acc) (const T *a, NUM v, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);
  if (!v || a->lo > a->hi) return;

  D *d = c->d;
  const NUM *ca = a->coef;
        NUM *cc = c->coef;
  ord_t new_hi = MIN3(a->hi,c->mo,d->trunc);
  ord_t new_lo = MIN(a->lo,c->lo);

  for (int i = d->ord2idx[new_lo ]; i < d->ord2idx[c->lo   ]; ++i) cc[i] = 0;
  for (int i = d->ord2idx[c->hi+1]; i < d->ord2idx[new_hi+1]; ++i) cc[i] = 0;
  for (int i = d->ord2idx[a->lo  ]; i < d->ord2idx[new_hi+1]; ++i) cc[i] += v * ca[i];

  c->lo = new_lo;
  c->hi = MAX(new_hi, c->hi);
  c->nz = mad_bit_trunc(mad_bit_add(c->nz,a->nz),c->hi);
}

/* TODO: check if this version is faster or not than the one above (use LINOP)
void
FUN(acc) (const T *a, NUM v, T *c)
{
  assert(a && c);
  ensure(a->d == c->d);
  if (!v || a->lo > a->hi) return;

  const T *t=0, *b=c;
  if (a->lo > b->lo) SWAP(a,b,t);

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->trunc);
  if (t) TPSA_LINOP(v*,+  );  // c->coef[i] = v*a->coef[i] +   c->coef[i];
  else   TPSA_LINOP(  ,+v*);  // c->coef[i] =   c->coef[i] + v*a->coef[i];
  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz), c->hi);
}
*/

void
FUN(add) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d);

  const T* t=0; 
  if (a->lo > b->lo) SWAP(a,b,t);

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->trunc);
  TPSA_LINOP( ,+);  // c->coef[i] = a->coef[i] + b->coef[i];
  c->lo = a->lo;    // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz), c->hi);
}

void
FUN(sub) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d);

  const T* t=0; 
  if (a->lo > b->lo) SWAP(a,b,t);

  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->trunc);
  if (t) TPSA_LINOP(-,+); // c->coef[i] = - a->coef[i] + b->coef[i];
  else   TPSA_LINOP( ,-); // c->coef[i] =   a->coef[i] - b->coef[i];
  c->lo = a->lo; // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz), c->hi);
}

void
FUN(mul) (const T *a, const T *b, T *r)
{
  assert(a && b && r);
  ensure(a->d == b->d && a->d == r->d);

  T *c = (a == r || b == r) ? r->d->PFX(t[0]) : r;

  D *d = a->d;
  c->lo = a->lo + b->lo;
  c->hi = MIN3(a->hi + b->hi, c->mo, d->trunc);

  // empty
  if (c->lo > c->hi) { FUN(clear)(c); goto ret; }

  // order 0
  NUM a0 = a->coef[0], b0 = b->coef[0];
  c->coef[0] = a0 * b0;
  c->nz = c->coef[0] != 0;
  if (c->hi == 0) {
    if (!c->coef[0]) c->lo = c->mo; // reset
    goto ret;
  }

  // order 1+
  idx_t max_ord1 = d->ord2idx[2];
  if (mad_bit_get(a->nz,1) && b0 && mad_bit_get(b->nz,1) && a0) {
    for (int i = 1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i] + b0*a->coef[i];
    c->nz = mad_bit_set(c->nz,1);
  }
  else if (mad_bit_get(a->nz,1) && b0) {
    for (int i = 1; i < max_ord1; ++i) c->coef[i] =                 b0*a->coef[i];
    c->nz = mad_bit_set(c->nz,1);
  }
  else if (mad_bit_get(b->nz,1) && a0) {
    for (int i = 1; i < max_ord1; ++i) c->coef[i] = a0*b->coef[i];
    c->nz = mad_bit_set(c->nz,1);
  }

  // order 2+
  if (c->hi >= 2) {
    if (a->lo > b->lo) {     //  a is the left-most one
      const T* t; SWAP(a,b,t);
      a0 = a->coef[0];
      b0 = b->coef[0];
    }

    ord_t c_hi = c->hi;              // needed by TPSA_LINOP
    TPSA_LINOP_ORD(b0 *, + a0 *, 2); // c->coef[i] = b0 * a->coef[i] + a0 * b->coef[i] ;
    for (int i = d->ord2idx[MAX(a->hi,b->hi)+1]; i < d->ord2idx[c_hi+1]; ++i)
      c->coef[i] = 0;

    if (a0) c->nz = mad_bit_trunc(mad_bit_add(c->nz,b->nz),c->hi);
    if (b0) c->nz = mad_bit_trunc(mad_bit_add(c->nz,a->nz),c->hi);

    #ifdef _OPENMP
    if (c->hi >= 12)
      hpoly_mul_par(a,b,c);
    else
    #endif
      hpoly_mul_ser(a,b,c);
  }

ret:
  assert(a != c && b != c);
  if (c != r) FUN(copy)(c,r);
}

void
FUN(div) (const T *a, const T *b, T *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && a->d == c->d);
  ensure(b->coef[0] != 0);

  if (b->hi == 0) { FUN(scl) (a,1/b->coef[0],c); return; }

  T *tmp = c->d->PFX(t[4]);  // t1-t3 used in inv
  FUN(inv) (b,1,tmp);
  FUN(mul) (a,tmp,c);
}

void
FUN(ipow) (const T *a, T *c, int n)
{
  assert(a && c);
  ensure(a->d == c->d);

  int inv = 0;

  if (n < 0) { n = -n; inv = 1; }

  T *t1 = c->d->PFX(t[1]);

  switch (n) {
    case 0: FUN(scalar) (c, 1);    break; // ok: no copy
    case 1: FUN(copy  ) (a, c);    break; // ok: 1 copy
    case 2: FUN(mul   ) (a,a, c);  break; // ok: 1 copy if a==c
    case 3: FUN(mul   ) (a,a, t1); FUN(mul)(t1,a,  c); break; // ok: 1 copy if a==c
    case 4: FUN(mul   ) (a,a, t1); FUN(mul)(t1,t1, c); break; // ok: no copy
    default: {
      T *t2 = c->d->PFX(t[2]);

      FUN(copy  )(a, t1);
      FUN(scalar)(c, 1 );

      for (;;) {
        if (n  & 1)   FUN(mul)(c ,t1, c ); // ok: 1 copy
        if (n /= 2) { FUN(mul)(t1,t1, t2); T *t=t2; t2=t1; t1=t; } // ok: no copy
        else break;
      }
    }
  }

  if (inv) FUN(inv)(c,1, c);
}

void
FUN(axpb) (NUM a, const T *x, NUM b, T *r)
{
  assert(x && r);
  ensure(x->d == r->d);
  FUN(scl)(x,a,r);
  if (b) FUN(set0)(r, 1,b);
}

void
FUN(axpbypc) (NUM c1, const T *a, NUM c2, const T *b, NUM c3, T *c)
{
  assert(a && b && c);
  ensure(a->d == b->d && b->d == c->d);

  if (a->lo > b->lo)  {
    const T* t; SWAP(a,b,t);
    NUM n;    SWAP(c1,c2,n);
  }
  ord_t   hi = MAX(a->hi,b->hi);
  ord_t c_hi = MIN3(hi, c->mo, c->d->trunc);  // TODO: optimise c_hi == 0 ?
  TPSA_LINOP(c1 *, + c2 *);  // c->coef[i] = c1 * a->coef[i] + c2 * b->coef[i];

  c->lo = a->lo;    // a->lo <= b->lo  (because of swap)
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz), c->hi);

  if (c3) FUN(set0) (c, 1,c3);
}

void
FUN(axypb) (NUM a, const T *x, const T *y, NUM b, T *r)
{
  assert(x && y && r);
  ensure(x->d == y->d && y->d == r->d);

  T *t1 = (x == r || y == r) ? r->d->PFX(t[1]) : r;
  FUN(mul)(x,y, t1);
  FUN(axpb)(a,t1, b, r);
}

void
FUN(axypbzpc) (NUM a, const T *x, const T *y, NUM b, const T *z, NUM c, T *r)
{
  assert(x && y && z && r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d);

  T *t1 = (x == r || y == r || z == r) ? r->d->PFX(t[1]) : r;
  FUN(mul)(x,y, t1);
  FUN(axpbypc)(a,t1, b,z, c, r);
}

void
FUN(axypbvwpc) (NUM a, const T *x, const T *y,
                NUM b, const T *v, const T *w, NUM c, T *r)
{
  assert(x && y && v && w && r);
  ensure(x->d == y->d && y->d == v->d && v->d == w->d && w->d == r->d);

  T *t1 = (x == r || y == r || v == r || w == r) ? r->d->PFX(t[1]) : r;
  T *t2 = (v == r || w == r || t1 == r) ? r->d->PFX(t[2]) : r;
  FUN(mul)(x,y, t1);
  FUN(mul)(v,w, t2);
  FUN(axpbypc)(a,t1, b,t2, c, r);
}

void
FUN(ax2pby2pcz2) (NUM a, const T *x, NUM b, const T *y, NUM c, const T *z, T *r)
{
  assert(x && y && z && r);
  ensure(x->d == y->d && y->d == z->d && z->d == r->d);

  T *t3 = (z == r) ? r->d->PFX(t[3]) : r;
  FUN(axypbvwpc)(a,x,x, b,y,y, 0, t3);
  FUN(axypbzpc)(c,z,z, 1,t3, 0, r);
}

void
FUN(poisson) (const T *a, const T *b, T *c, int n)
{
  // C = [A,B] (POISSON BRACKET, 2*n: No of PHASEVARS
  assert(a && b && c);
  ensure(a->d == b->d && b->d == c->d);

  T *is[4];
  for (int i = 0; i < 4; ++i)
    is[i] = FUN(new)(a, a->d->trunc);

  for (int i = 1; i <= n; ++i) {
    FUN(der)(a, is[0], 2*i - 1);
    FUN(der)(b, is[1], 2*i    );
    FUN(mul)(is[0],is[1],is[2]);
    FUN(add)(is[3],is[2],is[0]);
    FUN(copy)(is[0],is[3]);      // TODO: use swap ?

    FUN(der)(a, is[0], 2*i    );
    FUN(der)(b, is[1], 2*i - 1);
    FUN(mul)(is[0],is[1],is[2]);
    FUN(sub)(is[3],is[2],is[0]);
    FUN(copy)(is[0],is[3]);
  }

  FUN(copy)(is[3], c);
  for (int i = 0; i < 4; ++i)
    FUN(del)(is[i]);
}

// --- WITHOUT COMPLEX-BY-VALUE VERSION ---------------------------------------

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

#endif