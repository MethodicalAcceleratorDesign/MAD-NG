/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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
  ord_t o   = 0;
  idx_t i   = -1;
  log_t ok  = TRUE;
  idx_t *pi = t->d->ord2idx;

  ok &= t->mo <= t->d->mo && t->mo > 0;
  ok &= t->hi <= t->mo;
  ok &= t->lo <= t->hi || (t->lo == t->mo && t->hi == 0);
  ok &= t->lo ? !t->coef[0] : !!t->coef[0];

  if (!ok) goto ret;

  for (; o < t->lo; ++o) {
    ok &= !mad_bit_get(t->nz,o);
    if (!ok) goto ret;
  }
  for (; o <= t->hi; ++o) {
    if (!mad_bit_get(t->nz,o))
      for (i = pi[o]; i < pi[o+1]; ++i) {
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
  if (i_) *i_ = i < pi[o] ? pi[o] : i;
  return ok;
}

log_t
FUN(is_valid) (const T *t)
{
 assert(t);
 return FUN(check)(t, 0, 0);
}


void
FUN(debug) (const T *t, str_t name_, FILE *stream_)
{
  assert(t);
  if (!stream_) stream_ = stdout;

  fprintf(stream_, "{ %s: lo=%d hi=%d mo=%d",
          name_ ? name_ : "??", t->lo,t->hi,t->mo);

  ord_t  mo = MAX(0,MIN(t->mo,t->d->mo)); // avoid segfault if mo is corrupted
  idx_t *pi = t->d->ord2idx;
  for (idx_t i = 0; i < pi[mo+1]; ++i)
    fprintf(stream_," [%d]=" FMT, i, VAL(t->coef[i]));
  fprintf(stream_," }");

  char bnz[sizeof(t->nz)*CHAR_BIT+1] = {0};
  for (int b = 0; b <= mo; ++b)
    bnz[b] = '0' + !!mad_bit_get(t->nz,b);
  fprintf(stream_," nz=%s", bnz);

  ord_t o  = 0;
  idx_t i  = -1;
  log_t ok = FUN(check)(t, &o, &i);
  if (!ok) fprintf(stream_," ** (o=%d, i=%d)", o, i);
  fprintf(stream_,"\n");
}

// --- introspection ----------------------------------------------------------o

const D*
FUN(desc) (const T *t)
{
  assert(t);
  return t->d;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t);
  return t->d->ord2idx[t->mo+1];
}

ord_t
FUN(ord) (const T *t)
{
  assert(t);
  return t->mo;
}

ord_t
(FUN(ordv)) (const T *t1, const T *t2, ...)
{
  assert(t1 && t2);
  ord_t mo = MAX(t1->mo, t2->mo);
  va_list args;
  va_start(args,t2);
  while ((t2 = va_arg(args,T*)))
    if (t2->mo > mo) mo = t2->mo;
  va_end(args);
  return mo;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  if (!d) d = mad_desc_curr;
  ensure(d, "GTPSA descriptor not found");

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo, "GTPSA order exceeds descriptor maximum order");

  ssz_t nc = mad_desc_ordlen(d, mo);
  T *t = mad_malloc(sizeof(T) + nc * sizeof(NUM));
  t->d = d, t->mo = mo;
  return FUN(reset0)(t);
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t);
  return FUN(newd)(t->d, mo == mad_tpsa_same ? t->mo : mo);
}

void
FUN(del) (const T *t)
{
  mad_free((void*)t);
}

// --- clear, scalar ----------------------------------------------------------o

void
FUN(clear) (T *t)
{
  assert(t);
  FUN(reset0)(t);
}

void
FUN(scalar) (T *t, NUM v, idx_t iv)
{
  assert(t);
  if (!v && !iv) { FUN(reset0)(t); return; }

  t->lo = 0;
  t->coef[0] = v;

  if (iv > 0 && t->mo > 0) {
    ensure(iv <= t->d->nmv, "index exceeds GPTSA number of map variables");
    if (t->hi && mad_bit_get(t->nz,1)) {
      idx_t *pi = t->d->ord2idx;
      for (idx_t c = pi[1]; c < pi[2]; ++c) t->coef[c] = 0;
    }
    t->hi = 1;
    t->nz = 3;
    t->coef[iv] = 1;
  } else {
    t->hi = 0;
    t->nz = 1;
  }
}

// --- copy, convert ----------------------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r);
  if (t == r) return;
  const D *d = t->d;
  ensure(d == r->d, "incompatible GTPSAs descriptors");

  if (d->to < t->lo) { FUN(reset0)(r); return; }

  FUN(copy0)(t, r);

  idx_t *pi = d->ord2idx;
  for (idx_t i = pi[r->lo]; i < pi[r->hi+1]; ++i)
    r->coef[i] = t->coef[i];

  CHECK_VALIDITY(r);
}

void
FUN(convert) (const T *t, T *r, ssz_t n, idx_t t2r_[n])
{
  assert(t && r);
  if (t->d == r->d && !t2r_) { FUN(copy)(t,r); return; }
  ensure(t != r, "invalid destination (aliased source)");

  FUN(reset0)(r);

  ssz_t tn = t->d->nv, rn = r->d->nv;
  ord_t tm[tn], rm[rn];
  idx_t t2r[rn];

  idx_t i = 0;
  if (!t2r_)
    for (; i < MIN(rn,tn); ++i) t2r[i] = i; // identity
  else
    for (; i < MIN(rn, n); ++i)
      t2r[i] = t2r_[i] >= 0 && t2r_[i] < tn ? t2r_[i] : -1; // -1 -> discard var
  rn = i; // truncate

  idx_t *pi = t->d->ord2idx;
  for (idx_t ti = pi[t->lo]; ti < pi[t->hi+1]; ++ti) {
    if (!t->coef[ti]) continue;
    mad_desc_get_mono(t->d, tn, tm, ti);
    for (idx_t i = 0; i < rn; ++i)
      rm[i] = t2r[i] >= 0 ? tm[t2r[i]] : 0;
    if (mad_desc_mono_isvalid_m(r->d, rn, rm))
      FUN(setm)(r, rn, rm, 1, t->coef[ti]);
  }

  CHECK_VALIDITY(r);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, ssz_t n, ord_t m_[n], idx_t i)
{
  assert(t);
  return mad_desc_get_mono(t->d, n, m_, i);
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s);
  return mad_desc_get_idx_s(t->d, n, s);
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m);
  return mad_desc_get_idx_m(t->d, n, m);
}

idx_t
FUN(idxsm) (const T *t, ssz_t n, const idx_t m[n])
{
  assert(t && m);
  return mad_desc_get_idx_sm(t->d, n, m);
}

// --- accessors --------------------------------------------------------------o

NUM
FUN(get0) (const T *t)
{
  assert(t);
  return t->coef[0];
}

NUM
FUN(geti) (const T *t, idx_t i)
{
  assert(t);
  const D *d = t->d;
  ensure(i >= 0 && i < d->nc, "index order exceeds GPTSA maximum order");
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

void
FUN(getv) (const T *t, idx_t i, ssz_t n, NUM v[n])
{
  assert(t && v);
  const D *d = t->d;
  ensure(i >= 0 && i+n <= d->nc, "index order exceeds GPTSA maximum order");

  ord_t *ords = d->ords+i;
  const NUM *coef = t->coef+i;
  for (idx_t j = 0; j < n; ++j)
    v[j] = t->lo <= ords[j] && ords[j] <= t->hi ? coef[j] : 0;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{
  // --- mono is a string; represented as "[0-9]*"
  assert(t && s);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_s(d,n,s);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_m(d,n,m);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i,o)]
  assert(t && m);
  const D *d = t->d;
  idx_t i = mad_desc_get_idx_sm(d,n,m);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t);

  t->coef[0] = a*t->coef[0] + b;
  if (t->coef[0]) {
    idx_t *pi = t->d->ord2idx;
    for (idx_t c = pi[1]; c < pi[t->lo]; ++c) t->coef[c] = 0;
    t->nz = mad_bit_set(t->nz,0);
    t->lo = 0;
  }
  else {
    t->nz = mad_bit_clr(t->nz,0);
    int n = mad_bit_lowest(t->nz);
    t->lo = MIN(n,t->mo);
  }

  CHECK_VALIDITY(t);
}

void
FUN(seti) (T *t, idx_t i, NUM a, NUM b)
{
  assert(t);
  if (!i) { FUN(set0)(t,a,b); return; }

  const D *d = t->d;
  ensure(i > 0 && i < d->nc , "index order exceeds GPTSA maximum order");
  ensure(d->ords[i] <= t->mo, "index order exceeds GTPSA order");

  idx_t *pi = d->ord2idx;
  ord_t  o  = d->ords[i];

  NUM v = t->lo <= o && o <= t->hi ? a*t->coef[i]+b : b;

  if (!v) {
    t->coef[i] = v;
    idx_t j; // scan hpoly for non-zero
    for (j = pi[o]; j < pi[o+1] && !t->coef[j]; ++j);
    if (j == pi[o+1]) { // zero hpoly
      t->nz = mad_bit_clr(t->nz,o);
      int n = mad_bit_lowest(t->nz);
      t->lo = MIN(n,t->mo);
    }
    CHECK_VALIDITY(t);
    return;
  }

  if (t->lo > t->hi) {    // new TPSA, init ord o
    for (idx_t c = pi[o]; c < pi[o+1]; ++c) t->coef[c] = 0;
    t->lo = t->hi = o;
  }
  else if (o > t->hi) {   // extend right
    for (idx_t c = pi[t->hi+1]; c < pi[o+1]; ++c) t->coef[c] = 0;
    t->hi = o;
  }
  else if (o < t->lo) {   // extend left
    for (idx_t c = pi[o]; c < pi[t->lo]; ++c) t->coef[c] = 0;
    t->lo = o;
  }
  t->nz = mad_bit_set(t->nz,o);
  t->coef[i] = v;
  CHECK_VALIDITY(t);
}

void
FUN(setv) (T *t, idx_t i, ssz_t n, const NUM v[n])
{
  assert(t && v);
  const D *d = t->d;
  ensure(i >= 0 && i+n <= d->nc, "index order exceeds GPTSA maximum order");

  for (idx_t j = 0; j < n; j++) t->coef[j+i] = v[j];

  ord_t *ords = d->ords+i;
  if (t->lo > ords[0  ]) t->lo = ords[0  ];
  if (t->hi < ords[n-1]) t->hi = ords[n-1];

  idx_t *pi = d->ord2idx;
  for (idx_t o = ords[0]; o <= ords[n-1]; o++) {
    idx_t c; // scan hpoly for non-zero
    for (c = pi[o]; c < pi[o+1] && !t->coef[c]; ++c);
    t->nz = c == pi[o+1] ? mad_bit_clr(t->nz,o) : mad_bit_set(t->nz,o);
  }
  CHECK_VALIDITY(t);
}

void
FUN(sets) (T *t, ssz_t n, str_t s, NUM a, NUM b)
{
  assert(t && s);
  idx_t i = mad_desc_get_idx_s(t->d,n,s);
  FUN(seti)(t,i,a,b);
}

void
FUN(setm) (T *t, ssz_t n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m);
  idx_t i = mad_desc_get_idx_m(t->d,n,m);
  FUN(seti)(t,i,a,b);
}

void
FUN(setsm) (T *t, ssz_t n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m);
  idx_t i = mad_desc_get_idx_sm(t->d,n,m);
  FUN(seti)(t,i,a,b);
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

void FUN(scalar_r) (T *t, num_t v_re, num_t v_im, idx_t iv)
{ FUN(scalar)(t, CNUM(v), iv); }

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
