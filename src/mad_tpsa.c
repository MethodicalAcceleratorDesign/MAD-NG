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

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- debugging --------------------------------------------------------------o

void
FUN(debug) (const T *t)
{
  D *d = t->d;
  printf("{ nz=%d lo=%d hi=%d mo=%d | [0]=" FMT " ",
            t->nz,t->lo,t->hi,t->mo, VAL(t->coef[0]));

  ord_t hi = MIN3(t->hi, t->mo, t->d->trunc);
  int i = d->ord2idx[MAX(1,t->lo)]; // ord 0 already printed
  for (; i < d->ord2idx[hi+1]; ++i)
    if (t->coef[i])
      printf("[%d]=" FMT " ", i, VAL(t->coef[i]));
  printf(" }\n");
}

// --- introspection ----------------------------------------------------------o

D*
FUN(desc) (const T *t)
{
  assert(t);
  return t->d;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t);
  return mad_desc_tpsa_len(t->d, t->mo);
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
FUN(newd) (D *d, ord_t mo)
{
  assert(d);

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo, "GTPSA order exceeds descriptor maximum order");

  ssz_t nc = mad_desc_tpsa_len(d, mo);
  T *t = mad_malloc(sizeof(T) + nc * sizeof(NUM));

  t->d = d;
  t->lo = t->mo = mo;
  t->hi = t->nz = t->coef[0] = 0;  // coef[0] used without checking nz[0]
  return t;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t);
  if (mo == mad_tpsa_same) mo = t->mo;
  return FUN(newd)(t->d,mo);
}

void
FUN(copy) (const T *t, T *dst)
{
  assert(t && dst);
  ensure(t->d == dst->d, "incompatible GTPSAs descriptors");
  D *d = t->d;
  if (d->trunc < t->lo) { FUN(clear)(dst); return; }
  dst->hi = MIN3(t->hi, dst->mo, d->trunc);
  dst->lo = t->lo;
  dst->nz = mad_bit_trunc(t->nz, dst->hi);

  for (int i = d->ord2idx[dst->lo]; i < d->ord2idx[dst->hi+1]; ++i)
    dst->coef[i] = t->coef[i];
}

void
FUN(clear) (T *t)
{
  assert(t);
  t->hi = t->nz = t->coef[0] = 0;
  t->lo = t->mo;
}

void
FUN(scalar) (T *t, NUM v)
{
  assert(t);
  if (v) {
    t->coef[0] = v;
    t->nz = 1;
    t->lo = t->hi = 0;
  }
  else
    FUN(clear)(t);
}

#ifdef MAD_CTPSA_IMPL

void FUN(scalar_r) (T *t, num_t v_re, num_t v_im)
{ FUN(scalar)(t, CNUM(v)); }

#endif

void
FUN(del) (T *t)
{
  mad_free(t);
}

// --- indexing / monomials ---------------------------------------------------o

int
FUN(mono) (const T *t, int n, ord_t m_[n], idx_t i)
{
  assert(t);
  return mad_desc_get_mono(t->d, n, m_, i);
}

int
FUN(midx) (const T *t, int n, const ord_t m[n])
{
  assert(t && m);
  return mad_desc_get_idx(t->d, n, m);
}

int
FUN(midx_s) (const T *t, int n, str_t s)
{
  assert(t && s);
  return mad_desc_get_idx_s(t->d, n, s);
}

int
FUN(midx_sp) (const T *t, int n, const int m[n])
{
  assert(t && m);
  return mad_desc_get_idx_sp(t->d, n, m);
}

// --- accessors --------------------------------------------------------------o

NUM
FUN(get0) (const T *t)
{
  assert(t);
  return t->coef[0];
}

NUM
FUN(geti) (const T *t, int i)
{
  assert(t);
  D *d = t->d;
  ensure(i >= 0 && i < d->nc, "index out of bounds");
  ensure(d->ords[i] <= t->mo, "index order exceeds GPTSA maximum order");
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm) (const T *t, int n, const ord_t m[n])
{
  assert(t && m);
  D *d = t->d;
  idx_t i = mad_desc_get_idx(d,n,m);
  ensure(d->ords[i] <= t->mo, "monomial order exceeds GPTSA maximum order");
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm_s) (const T *t, int n, str_t s)
{
  // --- mono is a string; represented as "[0-9]*"
  assert(t && s);
  D *d = t->d;
  idx_t i = mad_desc_get_idx_s(d,n,s);
  ensure(d->ords[i] <= t->mo, "monomial order exceeds GPTSA maximum order");
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm_sp) (const T *t, int n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i o)]
  assert(t && m);
  D *d = t->d;
  idx_t i = mad_desc_get_idx_sp(d,n,m);
  ensure(d->ords[i] <= t->mo, "monomial order exceeds GPTSA maximum order");
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t);
  t->coef[0] = a*t->coef[0] + b;
  if (t->coef[0]) {
    t->nz = mad_bit_set(t->nz,0);
    for (int c = t->d->ord2idx[1]; c < t->d->ord2idx[t->lo]; ++c)
      t->coef[c] = 0;
    t->lo = 0;
  }
  else {
    int n = mad_bit_lowest(t->nz);
    t->nz = mad_bit_clr(t->nz,0);
    t->lo = MIN(n,t->mo);
  }
}

void
FUN(seti) (T *t, int i, NUM a, NUM b)
{
  assert(t);
  D *d = t->d;
  ensure(i >= 0 && i < d->nc, "index out of bounds");
  ensure(d->ords[i] <= t->mo, "index order exceeds GPTSA maximum order");

  if (i == 0) { FUN(set0)(t,a,b); return; }

  NUM v = a*FUN(geti)(t,i) + b;
  if (v == 0) {
    t->coef[i] = v;
    if (i == 0 && t->lo == 0) {
      int n = mad_bit_lowest(t->nz);
      t->nz = mad_bit_clr(t->nz,0);
      t->lo = MIN(n,t->mo);
    }
    return;
  }

  ord_t o = d->ords[i];
  t->nz = mad_bit_set(t->nz,o);
  if (t->lo > t->hi) {    // new TPSA, init ord o
    for (int c = d->ord2idx[o]; c < d->ord2idx[o+1]; ++c)
      t->coef[c] = 0;
    t->lo = t->hi = o;
  }
  else if (o > t->hi) {   // extend right
    for (int c = d->ord2idx[t->hi+1]; c < d->ord2idx[o+1]; ++c)
      t->coef[c] = 0;
    t->hi = o;
  }
  else if (o < t->lo) {   // extend left
    for (int c = d->ord2idx[o]; c < d->ord2idx[t->lo]; ++c)
      t->coef[c] = 0;
    t->lo = o;
  }
  t->coef[i] = v;
}

void
FUN(setm) (T *t, int n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m);
  assert(n <= t->d->nv);
  idx_t i = mad_desc_get_idx(t->d,n,m);
  FUN(seti)(t,i,a,b);
}

void
FUN(setm_s) (T *t, int n, str_t s, NUM a, NUM b)
{
  assert(t && s);
  idx_t i = mad_desc_get_idx_s(t->d,n,s);
  FUN(seti)(t,i,a,b);
}

void
FUN(setm_sp) (T *t, int n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m);
  idx_t i = mad_desc_get_idx_sp(t->d,n,m);
  FUN(seti)(t,i,a,b);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(get0_r) (const T *t, NUM *r)
{ assert(r); *r = FUN(get0)(t); }

void FUN(geti_r) (const T *t, int i, NUM *r)
{ assert(r); *r = FUN(geti)(t, i); }

void FUN(getm_r) (const T *t, int n, const ord_t m[n], NUM *r)
{ assert(r); *r = FUN(getm)(t, n, m); }

void FUN(getm_s_r) (const T *t, int n, str_t s, NUM *r)
{ assert(r); *r = FUN(getm_s)(t, n, s); }

void FUN(getm_sp_r) (const T *t, int n, const idx_t m[n], NUM *r)
{ assert(r); *r = FUN(getm_sp)(t, n, m); }

void FUN(set0_r) (T *t, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(set0)(t, CNUM(a), CNUM(b)); }

void FUN(seti_r) (T *t, int i, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(seti)(t, i, CNUM(a), CNUM(b)); }

void FUN(setm_r) (T *t, int n, const ord_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setm)(t, n, m, CNUM(a), CNUM(b)); }

void FUN(setm_s_r) (T *t, int n, str_t s, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setm_s)(t, n, s, CNUM(a), CNUM(b)); }

void FUN(setm_sp_r) (T *t, int n, const idx_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setm_sp)(t, n, m, CNUM(a), CNUM(b)); }

// --- end --------------------------------------------------------------------o

#endif // MAD_CTPSA_IMPL
