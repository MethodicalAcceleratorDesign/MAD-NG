/*
 o----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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

// --- DEBUGGING --------------------------------------------------------------

void
FUN(debug) (const T *t)
{
  D *d = t->desc;
  printf("{ nz=%d lo=%d hi=%d mo=%d | [0]=" FMT " ", t->nz, t->lo, t->hi, t->mo, VAL(t->coef[0]));
  ord_t hi = MIN3(t->hi, t->mo, t->desc->trunc);
  int i = d->hpoly_To_idx[MAX(1,t->lo)]; // ord 0 already printed
  for (; i < d->hpoly_To_idx[hi+1]; ++i)
    if (t->coef[i])
      printf("[%d]=" FMT " ", i, VAL(t->coef[i]));
  printf(" }\n");
}

// --- PUBLIC FUNCTIONS -------------------------------------------------------

// --- --- INTROSPECTION ------------------------------------------------------
D*
FUN(desc) (const T *t)
{
  assert(t);
  return t->desc;
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
  ord_t mo = t1->mo > t2->mo ? t1->mo : t2->mo;
  va_list args;
  va_start(args,t2);
  while ((t2 = va_arg(args,T*)))
    if (t2->mo > mo)
      mo = t2->mo;
  va_end(args);
  return mo;
}

// --- --- CTORS --------------------------------------------------------------

T*
FUN(newd) (D *d, ord_t mo)
{
  assert(d);

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo);

  T *t;
  if (d->PFX(stack_top) > 0)
    t = d->PFX(stack)[d->PFX(stack_top)--];
  else {
    t = malloc(sizeof(T) + d->nc * sizeof(NUM));
    assert(t);
    t->desc = d;
  }
  t->lo = t->mo = mo;
  t->hi = t->nz = t->coef[0] = 0;  // coef[0] used without checking NZ[0]
  t->tmp = 0;
  return t;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t);
  if (mo == mad_tpsa_same) mo = t->mo;
  return FUN(newd)(t->desc,mo);
}

T*
FUN(copy) (const T *src, T *dst)
{
  assert(src && dst);
  ensure(src->desc == dst->desc);
  D *d = src->desc;
  if (d->trunc < src->lo) {
    FUN(clear)(dst);
    return dst;
  }
  dst->hi = MIN3(src->hi, dst->mo, d->trunc);
  dst->lo = src->lo;
  dst->nz = mad_bit_trunc(src->nz, dst->hi);
  // dst->tmp = src->tmp;  // managed from outside

  for (int i = d->hpoly_To_idx[dst->lo]; i < d->hpoly_To_idx[dst->hi+1]; ++i)
    dst->coef[i] = src->coef[i];

  return dst;
}

void
FUN(clear) (T *t)
{
  assert(t);
  t->hi = t->nz = t->coef[0] = 0;
  t->lo = t->mo;
  // t->tmp = 0;  // managed from outside
}

void
FUN(del) (T *t)
{
  D *d = t->desc;
  if (d->PFX(stack_top) < mad_desc_stack)
    d->PFX(stack)[++d->PFX(stack_top)] = t;
  else
    free(t);
}

void
(FUN(delv))(T *t1, T *t2, ...)
{
  assert(t1 && t2);
  FUN(del)(t1);
  FUN(del)(t2);
  va_list args;
  va_start(args,t2);
  while((t2 = va_arg(args,T*)))
    FUN(del)(t2);
  va_end(args);
}

// --- --- INDEXING / MONOMIALS -----------------------------------------------

const ord_t*
FUN(mono) (const T *t, int i, int *n, ord_t *total_ord_)
{
  assert(t && n);
  D *d = t->desc;
  ensure(0 <= i && i < d->nc);
  *n = d->nv;
  if (total_ord_)
    *total_ord_ = d->ords[i];
  return d->To[i];
}

int
FUN(midx) (const T *t, int n, const ord_t m[n])
{
  assert(t && t->desc);
  assert(n <= t->desc->nv);
  return mad_desc_get_idx(t->desc, n, m);
}

int
FUN(midx_sp) (const T *t, int n, const int m[n])
{
  assert(t && t->desc);
  return mad_desc_get_idx_sp(t->desc, n, m);
}

// --- --- ACCESSORS ----------------------------------------------------------

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
  D *d = t->desc;
  ensure(i >= 0 && i < d->nc);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm) (const T *t, int n, const ord_t m[n])
{
  assert(t && m);
  D *d = t->desc;
  idx_t i = mad_desc_get_idx(d,n,m);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

NUM
FUN(getm_sp) (const T *t, int n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i o)]
  assert(t && m);
  D *d = t->desc;
  idx_t i = mad_desc_get_idx_sp(d,n,m);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
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
    FUN(clear) (t);
}

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t);
  t->coef[0] = a*t->coef[0] + b;
  if (t->coef[0]) {
    t->nz = mad_bit_set(t->nz,0);
    for (int c = t->desc->hpoly_To_idx[1]; c < t->desc->hpoly_To_idx[t->lo]; ++c)
      t->coef[c] = 0;
    t->lo = 0;
  }
  else {
    t->nz = mad_bit_clr(t->nz,0);
    t->lo = MIN(mad_bit_lowest(t->nz),t->mo);
  }
}

void
FUN(seti) (T *t, int i, NUM a, NUM b)
{
  assert(t);
  D *d = t->desc;
  ensure(i >= 0 && i < d->nc && d->ords[i] <= t->mo);

  if (i == 0) { FUN(set0)(t,a,b); return; }

  NUM v = a*FUN(geti)(t,i) + b;
  if (v == 0) {
    t->coef[i] = v;
    if (i == 0 && t->lo == 0) {
      t->nz = mad_bit_clr(t->nz,0);
      t->lo = MIN(mad_bit_lowest(t->nz),t->mo);
    }
    return;
  }

  ord_t o = d->ords[i];
  t->nz = mad_bit_set(t->nz,o);
  if (t->lo > t->hi) {    // new TPSA, init ord o
    for (int c = d->hpoly_To_idx[o]; c < d->hpoly_To_idx[o+1]; ++c)
      t->coef[c] = 0;
    t->lo = t->hi = o;
  }
  else if (o > t->hi) {   // extend right
    for (int c = d->hpoly_To_idx[t->hi+1]; c < d->hpoly_To_idx[o+1]; ++c)
      t->coef[c] = 0;
    t->hi = o;
  }
  else if (o < t->lo) {   // extend left
    for (int c = d->hpoly_To_idx[o]; c < d->hpoly_To_idx[t->lo]; ++c)
      t->coef[c] = 0;
    t->lo = o;
  }
  t->coef[i] = v;
}

void
FUN(setm) (T *t, int n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m);
  assert(n <= t->desc->nv);
  idx_t i = mad_desc_get_idx(t->desc,n,m);
  FUN(seti)(t,i,a,b);
}

void
FUN(setm_sp) (T *t, int n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m);
  idx_t i = mad_desc_get_idx_sp(t->desc,n,m);
  FUN(seti)(t,i,a,b);
}

// --- --- TRANSFORMATION -----------------------------------------------------

T*
FUN(map) (const T *a, T *c, NUM (*f)(NUM v, int i_))
{
  assert(a && c);
  ensure(a->desc == c->desc);
  // TODO: use on the whole range, not just [lo,hi]

  D *d = a->desc;
  if (d->trunc < a->lo) { FUN(clear)(c); return c; }

  c->hi = MIN3(a->hi, c->mo, d->trunc);
  c->lo = a->lo;
  c->nz = mad_bit_trunc(a->nz, c->hi);

  for (int i = d->hpoly_To_idx[c->lo]; i < d->hpoly_To_idx[c->hi+1]; ++i)
    c->coef[i] = f(a->coef[i], i);

  return c;
}

T*
FUN(map2) (const T *a, const T *b, T *c, NUM (*f)(NUM va, NUM vb, int i_))
{
  assert(a && b && c);
  ensure(a->desc == b->desc && a->desc == c->desc);

  idx_t *pi = a->desc->hpoly_To_idx;
  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }
  ord_t c_hi = MIN3(MAX(a->hi,b->hi),c->mo,c->desc->trunc),
        c_lo = a->lo;

  NUM va, vb;
  for (int i = pi[c_lo]; i < pi[c_hi+1]; ++i) {
    int curr_strict_mode = mad_tpsa_strict;
    mad_tpsa_strict = 0;
    va = FUN(geti)(a,i);
    vb = FUN(geti)(b,i);
    mad_tpsa_strict = curr_strict_mode;
    c->coef[i] = f(va,vb,i);
  }
  c->lo = c_lo;
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz),c->hi);

  return c;
}
