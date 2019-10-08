/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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

#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <limits.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- debugging --------------------------------------------------------------o

/*
GTPSA are dense within [lo,hi], nz[ord] is zero outside, coef[0] is always set
ord in [lo,hi]: nz[ord] == 0 <=> all coef[ord] == 0
ord in [lo,hi]: nz[ord] == 1 <=> any coef[ord] != 0 (or none)
                nz[0  ] == 0 <=>     coef[0  ] == 0 && lo >= 1
                nz[0  ] == 1 <=>     coef[0  ] != 0 && lo == 0
GTPSA just initialized have lo=mo, hi=0, nz=0, coef[0]=0 (see reset)
*/

static inline log_t
FUN(check) (const T *t, ord_t *o_, idx_t *i_)
{
  assert(t); assert(t->d);
  ord_t o  = 0;
  idx_t i  = -1;
  log_t ok = TRUE;

  ok &= t->mo <= t->d->mo;
  ok &= t->hi <= t->mo;
  ok &= t->lo <= t->hi || (t->lo == t->mo && t->hi == 0);
  ok &= t->lo <= t->hi && !t->lo ? !!t->coef[0] : !t->coef[0];

  if (!ok) goto ret; // i == -1

  for (; o < t->lo; ++o) {
    ok &= !mad_bit_get(t->nz,o);
    if (!ok) goto ret;
  }

  idx_t *o2i = t->d->ord2idx;

  for (; o <= t->hi; ++o) {
    if (!mad_bit_get(t->nz,o))
      for (i = o2i[o]; i < o2i[o+1]; ++i) {
        ok &= !t->coef[i];
        if (!ok) goto ret;
      }
  }

  for (; o <= t->mo; ++o) {
    ok &= !mad_bit_get(t->nz,o);
    if (!ok) goto ret;
  }

ret:
  if (o_) *o_ = o;
  if (i_) *i_ = i;

  return ok;
}

log_t
FUN(is_valid) (const T *t)
{
 DBGFUN(->);
 log_t ret = FUN(check)(t,0,0);
 DBGFUN(<-); return ret;
}

void
FUN(debug) (const T *t, str_t name_, int line, FILE *stream_)
{
  assert(t); assert(t->d); DBGFUN(->);
  const D* d = t->d;
  if (!stream_) stream_ = stdout;

  fprintf(stream_, "{ %s:%d: lo=%d hi=%d mo=%d id=%d",
          name_ ? name_ : "??", line, t->lo,t->hi,t->mo,d->id); fflush(stream_);

  ord_t   mo = MIN(t->mo,d->mo); // avoid segfault if mo is corrupted
  idx_t *o2i = d->ord2idx;
  for (idx_t i = 0; i < o2i[mo+1]; ++i)
    fprintf(stream_," [%d]=" FMT, i, VAL(t->coef[i]));
  fprintf(stream_," }"); fflush(stream_);

  char bnz[sizeof(t->nz)*CHAR_BIT+1] = {0};
  for (ord_t b = 0; b <= mo; ++b)
    bnz[b] = '0' + !!mad_bit_get(t->nz,b);
  fprintf(stream_," nz=%s", bnz); fflush(stream_);

  ord_t o  = 0;
  idx_t i  = -1;
  log_t ok = FUN(check)(t, &o, &i);
  if (!ok) fprintf(stream_," ** (o=%d, i=%d)\n", o, i);
  else     fprintf(stream_,"\n");
  DBGFUN(<-);
}

// --- introspection ----------------------------------------------------------o

const D*
FUN(desc) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D* ret = t->d;
  DBGFUN(<-); return ret;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ssz_t ret = t->d->ord2idx[t->mo+1];
  DBGFUN(<-); return ret;
}

ord_t
FUN(ord) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ord_t ret = t->mo;
  DBGFUN(<-); return ret;
}

ord_t
(FUN(ordv)) (const T *t, ...)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ord_t mo = t->mo;
  va_list args;
  va_start(args,t);
  while ((t = va_arg(args,T*)))
    if (t->mo > mo) mo = t->mo;
  va_end(args);
  DBGFUN(<-); return mo;
}

ord_t
FUN(ordn) (ssz_t n, const T *t[])
{
  assert(t); DBGFUN(->);
  ord_t mo = 0;
  for (idx_t i = 0; i < n; i++)
    if (t[i]->mo > mo) mo = t[i]->mo;
  DBGFUN(<-); return mo;
}

// --- init (unsafe) ----------------------------------------------------------o

T*
FUN(init) (T *t, const D *d, ord_t lo, ord_t hi, ord_t mo)
{
  assert(t); DBGFUN(->);
  if (!d) { assert(mad_desc_curr); d = mad_desc_curr; }
  ensure(d, "GTPSA descriptor not found");

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo,
              "GTPSA order %d exceeds descriptor maximum order %d", mo, d->mo);

  t->d = d, t->lo = MIN(lo,mo), t->hi = MIN(hi,mo), t->mo = mo;

  DBGTPSA(t); DBGFUN(<-); return t;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  DBGFUN(->);
  if (!d) d = mad_desc_curr;
  ensure(d, "GTPSA descriptor not found (no current one?)");

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo, "GTPSA order exceeds descriptor maximum order");

  ssz_t nc = mad_desc_ordlen(d, mo);
  T *t = mad_malloc(sizeof(T) + nc * sizeof(NUM));
  t->d = d, t->mo = mo;
  FUN(reset0)(t);

  DBGTPSA(t); DBGFUN(<-); return t;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  T *ret = FUN(newd)(t->d, mo == mad_tpsa_same ? t->mo : mo);
  DBGFUN(<-); return ret;
}

void
FUN(del) (const T *t)
{
  DBGFUN(->);
  if (t) { DBGTPSA(t); mad_free((void*)t); }
  DBGFUN(<-);
}

// --- clear, scalar ----------------------------------------------------------o

void
FUN(clear) (T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  FUN(reset0)(t);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(scalar) (T *t, NUM v, idx_t iv, NUM scl)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!v && !(iv && t->mo)) { FUN(reset0)(t); DBGTPSA(t); DBGFUN(<-); return; }

  t->coef[0] = v;

  if (iv && t->mo) {
    const D *d = t->d;
    ensure(0 <= iv && iv <= d->nv,
           "index %d exceeds GPTSA number of variables %d", iv, d->nv);

    if (t->hi && mad_bit_get(t->nz,1)) {
      idx_t *o2i = d->ord2idx; // clear first order
      for (idx_t i = o2i[1]; i < o2i[2]; ++i) t->coef[i] = 0;
    }
    t->hi = 1;
    t->lo = v ? 0 : 1;
    t->nz = v ? 3 : 2;
    t->coef[iv] = scl ? scl : 1;
  } else {
    t->hi = 0;
    t->lo = v ? 0 : t->mo;
    t->nz = v ? 1 : 0;
  }

  DBGTPSA(t); DBGFUN(<-);
}

// --- copy, convert ----------------------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  if (t == r) return;
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);
  const D *d = t->d;

  if (d->to < t->lo) { FUN(reset0)(r); DBGTPSA(r); DBGFUN(<-); return; }

  FUN(copy0)(t, r);

  idx_t *o2i = d->ord2idx;
  for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i)
    r->coef[i] = t->coef[i];

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(getord) (const T *t, T *r, ord_t ord)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);
  const D *d = t->d;

  if (ord < t->lo || ord > MIN(t->hi, r->mo)) {
    FUN(reset0)(r); DBGTPSA(r); DBGFUN(<-); return;
  }

  r->lo = r->hi = ord;
  r->nz = mad_bit_lcut(t->nz,r->lo);
  r->nz = mad_bit_hcut(r->nz,r->hi);

  if (r != t) {
    idx_t *o2i = d->ord2idx;
    for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i)
      r->coef[i] = t->coef[i];
  }
  if (r->lo) r->coef[0] = 0;

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(cutord) (const T *t, T *r, int ord)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);
  const D *d = t->d;

  if (ord < 0) {
    ord_t o = MIN3(-ord, t->mo, d->to);
    if (o >= t->hi) { FUN(reset0)(r); DBGTPSA(r); DBGFUN(<-); return; }
    FUN(copy0)(t, r);
    r->lo = MIN(MAX(t->lo, o+1), r->mo);
    r->nz = mad_bit_lcut(t->nz,r->lo);
  } else {
    ord_t o = MIN3(ord, t->mo, d->to);
    if (o <= t->lo) { FUN(reset0)(r); DBGTPSA(r); DBGFUN(<-); return; }
    FUN(copy0)(t, r);
    r->hi = MIN3(t->hi, o-1, r->mo);
    r->nz = mad_bit_hcut(t->nz,r->hi);
  }

  if (r != t) {
    idx_t *o2i = d->ord2idx;
    for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i)
      r->coef[i] = t->coef[i];
  }
  if (r->lo) r->coef[0] = 0;

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(convert) (const T *t, T *r_, ssz_t n, idx_t t2r_[n])
{
  assert(t && r_); DBGFUN(->); DBGTPSA(t); DBGTPSA(r_);

  if (t->d == r_->d && !t2r_) { FUN(copy)(t,r_); DBGFUN(<-); return; }

  T *r = t == r_ ? GET_TMPX(r_) : r_;

  FUN(reset0)(r);

  ssz_t rn = r->d->nv, tn = t->d->nv;
  ord_t rm[rn], tm[tn];
  idx_t t2r[rn]; // rm[i] = tm[t2r[i]], i=0..rn-1

  idx_t i = 0;
  if (!t2r_)
    for (; i < MIN(rn,tn); ++i) t2r[i] = i; // identity
  else
    for (; i < MIN(rn, n); ++i)
      t2r[i] = t2r_[i] > 0 && t2r_[i] <= tn ? t2r_[i]-1 : -1; // -1 -> discard var
  rn = i; // truncate

  idx_t *o2i = t->d->ord2idx;
  for (idx_t ti = o2i[t->lo]; ti < o2i[t->hi+1]; ++ti) {
    if (!t->coef[ti]) goto skip;
    mad_desc_get_mono(t->d, tn, tm, ti);
    for (idx_t i = 0; i < rn; ++i) {
      if (tm[i] && t2r[i] < 0) goto skip;   // discard var
      rm[i] = tm[t2r[i]];                   // translate/permute
    }
    idx_t ir = mad_desc_get_idx_m(r->d, rn, rm);
    if (ir >= 0) FUN(seti)(r, ir, 0, t->coef[ti]);
  skip: ;
  }

  if (r != r_) { FUN(copy)(r,r_); REL_TMPX(r); }

  DBGTPSA(r_); DBGFUN(<-);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, ssz_t n, ord_t m_[n], idx_t i)
{
  assert(t); DBGTPSA(t);
  ord_t ret = mad_desc_get_mono(t->d, n, m_, i);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s); DBGTPSA(t);
  idx_t ret = mad_desc_get_idx_s(t->d, n, s);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGTPSA(t);
  idx_t ret = mad_desc_get_idx_m(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxsm) (const T *t, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t ret = mad_desc_get_idx_sm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(cycle) (const T *t, ssz_t n, ord_t m_[n], idx_t i, num_t *v_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  ensure(0 <= i && i < d->nc, "index out of bounds");

  idx_t *o2i = d->ord2idx;
  idx_t  ni = o2i[t->hi+1];
  for (i = MAX(i,o2i[t->lo]); i < ni && !t->coef[i]; ++i) ;

  if (i >= ni) { DBGFUN(<-); return -1; }

  if (m_) {
    ensure(0 <= n && n <= d->nv, "invalid monomial length");
    mad_mono_copy(n, d->To[i], m_);
  }
  if (v_) *v_ = t->coef[i];

  idx_t ret = ++i == ni ? -1 : i;
  DBGFUN(<-); return ret;
}

// --- accessors --------------------------------------------------------------o

NUM
FUN(get0) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  NUM ret = t->coef[0];
  DBGFUN(<-); return ret;
}

NUM
FUN(geti) (const T *t, idx_t i)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  ensure(i >= 0 && i < d->nc, "index order exceeds GPTSA maximum order");
  NUM ret = t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{
  // --- mono is a string; represented as "[0-9A-Z]*"
  assert(t && s); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_s(d,n,s);
  NUM ret = t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_m(d,n,m);
  NUM ret = t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i,o)]
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_sm(d,n,m);
  NUM ret = t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  t->coef[0] = a*t->coef[0]+b;
  if (t->coef[0]) {
    idx_t *o2i = t->d->ord2idx;
    for (idx_t c = o2i[1]; c < o2i[t->lo]; ++c) t->coef[c] = 0;
    t->nz = mad_bit_set(t->nz,0);
    t->lo = 0;
  }
  else {
    t->nz = mad_bit_clr(t->nz,0);
    t->lo = MIN(mad_bit_lowest(t->nz),t->mo);
    t->hi = t->nz ? MIN(mad_bit_highest(t->nz),t->mo) : 0;
  }

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(seti) (T *t, idx_t i, NUM a, NUM b)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!i) { FUN(set0)(t,a,b); DBGTPSA(t); DBGFUN(<-); return; }

  const D *d = t->d;
  ensure(i > 0 && i <  d->nc, "index order exceeds GPTSA maximum order");
  ensure(d->ords[i] <= t->mo, "index order exceeds GTPSA order");

  idx_t *o2i = d->ord2idx;
  ord_t  o   = d->ords[i];

  NUM v = t->lo <= o && o <= t->hi ? a*t->coef[i]+b : b;

  if (!v) {
    if (!mad_bit_get(t->nz,o)) { DBGTPSA(t); DBGFUN(<-); return; } // already zeroed
    t->coef[i] = 0;
    idx_t j; // scan hpoly for non-zero
    for (j = o2i[o]; j < o2i[o+1] && !t->coef[j]; ++j) ;
    if (j == o2i[o+1]) { // zero hpoly
      t->nz = mad_bit_clr(t->nz,o);
      t->lo = MIN(mad_bit_lowest(t->nz),t->mo);
      t->hi = t->nz ? MIN(mad_bit_highest(t->nz),t->mo) : 0;
      if (t->lo) t->coef[0] = 0;
    }
    DBGTPSA(t); DBGFUN(<-); return;
  }

  if (t->lo > t->hi) {    // new TPSA, init ord o
    for (idx_t c = o2i[o]; c < o2i[o+1]; ++c) t->coef[c] = 0;
    t->lo = t->hi = o;
  }
  else if (o > t->hi) {   // extend right
    for (idx_t c = o2i[t->hi+1]; c < o2i[o+1]; ++c) t->coef[c] = 0;
    t->hi = o;
  }
  else if (o < t->lo) {   // extend left
    for (idx_t c = o2i[o]; c < o2i[t->lo]; ++c) t->coef[c] = 0;
    t->lo = o;
  }
  t->nz = mad_bit_set(t->nz,o);
  t->coef[i] = v;

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(sets) (T *t, ssz_t n, str_t s, NUM a, NUM b)
{
  assert(t && s); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_get_idx_s(t->d,n,s);
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setm) (T *t, ssz_t n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_get_idx_m(t->d,n,m);
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setsm) (T *t, ssz_t n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_get_idx_sm(t->d,n,m);
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

// --- vector-based accessors -------------------------------------------------o

void
FUN(getv) (const T *t, idx_t i, ssz_t n, NUM v[n])
{
  assert(t && v); DBGFUN(->); DBGTPSA(t);

  const D *d = t->d;
  ensure(i >= 0 && i+n <= d->nc, "index order exceeds GPTSA maximum order");

  ord_t *ords = d->ords+i;
  const NUM *coef = t->coef+i;
  for (idx_t j = 0; j < n; ++j)
    v[j] = t->lo <= ords[j] && ords[j] <= t->hi ? coef[j] : 0;

  DBGFUN(<-);
}

void
FUN(setv) (T *t, idx_t i, ssz_t n, const NUM v[n])
{
  assert(t && v); DBGFUN(->); DBGTPSA(t);

  const D *d = t->d;
  ensure(i >= 0 && i+n  <= d->nc, "index order exceeds GPTSA maximum order");
  ensure(d->ords[i+n-1] <= t->mo, "index order exceeds GTPSA order");

  for (idx_t j = 0; j < n; j++) t->coef[j+i] = v[j];

  ord_t *ords = d->ords+i;
  if (t->lo > ords[0  ]) t->lo = ords[0  ];
  if (t->hi < ords[n-1]) t->hi = ords[n-1];

  idx_t *o2i = d->ord2idx;
  for (idx_t o = ords[0]; o <= ords[n-1]; o++) {
    idx_t c; // scan hpoly for non-zero
    for (c = o2i[o]; c < o2i[o+1] && !t->coef[c]; ++c) ;
    t->nz = c == o2i[o+1] ? mad_bit_clr(t->nz,o) : mad_bit_set(t->nz,o);
  }
  t->lo = MIN(mad_bit_lowest(t->nz),t->mo);
  t->hi = t->nz ? MIN(mad_bit_highest(t->nz),t->mo) : 0;
  if (t->lo) t->coef[0] = 0;

  DBGTPSA(t); DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(get0_r) (const T *t, NUM *r)
{ assert(r); *r = FUN(get0)(t); }

void FUN(geti_r) (const T *t, idx_t i, NUM *r)
{ assert(r); *r = FUN(geti)(t, i); }

void FUN(gets_r) (const T *t, ssz_t n, str_t s, NUM *r)
{ assert(r); *r = FUN(gets)(t, n, s); }

void FUN(getm_r) (const T *t, ssz_t n, const ord_t m[n], NUM *r)
{ assert(r); *r = FUN(getm)(t, n, m); }

void FUN(getsm_r) (const T *t, ssz_t n, const idx_t m[n], NUM *r)
{ assert(r); *r = FUN(getsm)(t, n, m); }

void FUN(set0_r) (T *t, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(set0)(t, CNUM(a), CNUM(b)); }

void FUN(seti_r) (T *t, idx_t i, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(seti)(t, i, CNUM(a), CNUM(b)); }

void FUN(sets_r) (T *t, ssz_t n, str_t s, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(sets)(t, n, s, CNUM(a), CNUM(b)); }

void FUN(setm_r) (T *t, ssz_t n, const ord_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setm)(t, n, m, CNUM(a), CNUM(b)); }

void FUN(setsm_r) (T *t, ssz_t n, const idx_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setsm)(t, n, m, CNUM(a), CNUM(b)); }

void FUN(scalar_r) (T *t, num_t v_re, num_t v_im, idx_t iv, num_t scl_re, num_t scl_im)
{ FUN(scalar)(t, CNUM(v), iv, CNUM(scl)); }

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
