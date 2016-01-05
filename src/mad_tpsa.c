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
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"

#include "mad_bit.h"
#include "mad_mono.h"
#include "mad_tpsa.h"

#undef  ensure
#define ensure(test) mad_ensure(test, MKSTR(test))

#include "mad_tpsa_desc.tc"

// #define TRACE

#define T struct tpsa
#define D struct tpsa_desc

const ord_t mad_tpsa_default   = -1;
const ord_t mad_tpsa_same      = -2;
      int   mad_tpsa_strict    =  0;

struct tpsa { // warning: must be kept identical to LuaJIT definition (cmad.lua)
  D      *desc;
  ord_t   lo, hi, mo; // lowest/highest used ord, trunc ord
  bit_t   nz;
  int     tmp;
  num_t   coef[];
};

// --- DEBUGGING --------------------------------------------------------------

void
mad_tpsa_debug(const T *t)
{
  D *d = t->desc;
  printf("{ nz=%d lo=%d hi=%d mo=%d | [0]=%g ", t->nz, t->lo, t->hi, t->mo, t->coef[0]);
  ord_t hi = min_ord(t->hi, t->mo, t->desc->trunc);
  int i = d->hpoly_To_idx[MAX(1,t->lo)]; // ord 0 already printed
  for (; i < d->hpoly_To_idx[hi+1]; ++i)
    if (t->coef[i])
      printf("[%d]=%g ", i, t->coef[i]);
  printf(" }\n");
}

// --- PUBLIC FUNCTIONS -------------------------------------------------------

// --- --- INTROSPECTION ------------------------------------------------------
D*
mad_tpsa_desc(const T *t)
{
  assert(t);
  return t->desc;
}

ord_t
mad_tpsa_gtrunc(D *d, ord_t to)
{
  assert(d);
  ord_t orig = d->trunc;
  if (to == mad_tpsa_same)
    return orig;

  if (to == mad_tpsa_default)
    to = d->mo;
  else
    ensure(to <= d->mo);
  d->trunc = to;
  return orig;
}

ord_t
mad_tpsa_ord(const T *t)
{
  assert(t);
  return t->mo;
}

ord_t
(mad_tpsa_ordv)(const T *t1, const T *t2, ...)
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

int
mad_tpsa_maxsize(const D *d)
{
  assert(d);
  return d->nc;
}

ord_t
mad_tpsa_maxord(const D *d)
{
  assert(d);
  return d->mo;
}

// --- --- CTORS --------------------------------------------------------------

T*
mad_tpsa_newd(D *d, ord_t mo)
{
  assert(d);

  if (mo == mad_tpsa_default)
    mo = d->mo;
  else
    ensure(mo <= d->mo);

  T *t;
  if (d->stack_top > 0)
    t = d->stack[d->stack_top--];
  else {
    t = malloc(sizeof(T) + d->nc * sizeof(num_t));
    assert(t);
    t->desc = d;
  }
  t->lo = t->mo = mo;
  t->hi = t->nz = t->coef[0] = 0;  // coef[0] used without checking NZ[0]
  t->tmp = 0;
  return t;
}

T*
mad_tpsa_new(const T *t, ord_t mo)
{
  assert(t);
  if (mo == mad_tpsa_same)
    mo = t->mo;
  return mad_tpsa_newd(t->desc,mo);
}

T*
mad_tpsa_copy(const T *src, T *dst)
{
  assert(src && dst);
  ensure(src->desc == dst->desc);
  D *d = src->desc;
  if (d->trunc < src->lo) {
    mad_tpsa_clear(dst);
    return dst;
  }
  dst->hi = min_ord(src->hi, dst->mo, d->trunc);
  dst->lo = src->lo;
  dst->nz = mad_bit_trunc(src->nz, dst->hi);
  // dst->tmp = src->tmp;  // managed from outside

  for (int i = d->hpoly_To_idx[dst->lo]; i < d->hpoly_To_idx[dst->hi+1]; ++i)
    dst->coef[i] = src->coef[i];
#ifdef TRACE
  printf("Copied from %p to %p\n", (void*)src, (void*)dst);
#endif
  return dst;
}

void
mad_tpsa_clear(T *t)
{
  assert(t);
  t->hi = t->nz = t->coef[0] = 0;
  t->lo = t->mo;
  // t->tmp = 0;  // managed from outside
}

void
mad_tpsa_del(T *t)
{
#ifdef TRACE
  printf("tpsa del %p\n", (void*)t);
#endif
  D *d = t->desc;
  if (d->stack_top < DESC_STACK_SIZE)
    d->stack[++d->stack_top] = t;
  else
    free(t);
}

void
(mad_tpsa_delv)(T *t1, T *t2, ...)
{
  assert(t1 && t2);
  mad_tpsa_del(t1);
  mad_tpsa_del(t2);
  va_list args;
  va_start(args,t2);
  while((t2 = va_arg(args,T*)))
    mad_tpsa_del(t2);
  va_end(args);
}

// --- --- INDEXING / MONOMIALS -----------------------------------------------

const ord_t*
mad_tpsa_mono(const T *t, int i, int *n, ord_t *total_ord_)
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
mad_tpsa_midx(const T *t, int n, const ord_t m[n])
{
  assert(t && t->desc);
  assert(n <= t->desc->nv);
  return desc_get_idx(t->desc, n, m);
}

int
mad_tpsa_midx_sp(const T *t, int n, const int m[n])
{
  assert(t && t->desc);
  return desc_get_idx_sp(t->desc, n, m);
}


// --- --- ACCESSORS ----------------------------------------------------------

num_t
mad_tpsa_get0(const T *t)
{
  assert(t);
  return t->coef[0];
}

num_t
mad_tpsa_geti(const T *t, int i)
{
  assert(t);
  D *d = t->desc;
  ensure(i >= 0 && i < d->nc);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

num_t
mad_tpsa_getm(const T *t, int n, const ord_t m[n])
{
  assert(t && m);
  D *d = t->desc;
  idx_t i = desc_get_idx(d,n,m);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

num_t
mad_tpsa_getm_sp(const T *t, int n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i o)]
  assert(t && m);
  D *d = t->desc;
  idx_t i = desc_get_idx_sp(d,n,m);
  if (mad_tpsa_strict)
    ensure(d->ords[i] <= t->mo);
  return t->lo <= d->ords[i] && d->ords[i] <= t->hi ? t->coef[i] : 0;
}

void
mad_tpsa_scalar(T *t, num_t v)
{
  assert(t);
  if (v) {
    t->coef[0] = v;
    t->nz = 1;
    t->lo = t->hi = 0;
  }
  else
    mad_tpsa_clear(t);
}

void
mad_tpsa_set0(T *t, num_t a, num_t b)
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
    t->lo = min_ord2(mad_bit_lowest(t->nz),t->mo);
  }
}

void
mad_tpsa_seti(T *t, int i, num_t a, num_t b)
{
#ifdef TRACE
  printf("tpsa_seti for %p i=%d a=%lf b=%lf\n", (void*)t, i, a,b);
#endif
  assert(t);
  D *d = t->desc;
  ensure(i >= 0 && i < d->nc && d->ords[i] <= t->mo);

  if (i == 0) { mad_tpsa_set0(t,a,b); return; }

  num_t v = a*mad_tpsa_geti(t,i) + b;
  if (v == 0) {
    t->coef[i] = v;
    if (i == 0 && t->lo == 0) {
      t->nz = mad_bit_clr(t->nz,0);
      t->lo = min_ord2(mad_bit_lowest(t->nz),t->mo);
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
mad_tpsa_setm(T *t, int n, const ord_t m[n], num_t a, num_t b)
{
  assert(t && m);
  assert(n <= t->desc->nv);
#ifdef TRACE
  printf("set_mono: "); mad_mono_print(n, m); printf("\n");
#endif
  idx_t i = desc_get_idx(t->desc,n,m);
  mad_tpsa_seti(t,i,a,b);
}

void
mad_tpsa_setm_sp(T *t, int n, const idx_t m[n], num_t a, num_t b)
{
  assert(t && m);
#ifdef TRACE
  printf("set_mono_sp: [ ");
  for (int i=0; i < n; i += 2)
    printf("%d %d  ", m[i], m[i+1]);
  printf("]\n");
#endif
  idx_t i = desc_get_idx_sp(t->desc,n,m);
  mad_tpsa_seti(t,i,a,b);
}

// --- --- TRANSFORMATION -----------------------------------------------------

T*
mad_tpsa_map(const T *a, T *c, num_t (*f)(num_t v, int i_))
{
  assert(a && c);
  ensure(a->desc == c->desc);
  // TODO: use on the whole range, not just [lo,hi]

  D *d = a->desc;
  if (d->trunc < a->lo) { mad_tpsa_clear(c); return c; }

  c->hi = min_ord(a->hi, c->mo, d->trunc);
  c->lo = a->lo;
  c->nz = mad_bit_trunc(a->nz, c->hi);

  for (int i = d->hpoly_To_idx[c->lo]; i < d->hpoly_To_idx[c->hi+1]; ++i)
    c->coef[i] = f(a->coef[i], i);

  return c;
}

T*
mad_tpsa_map2(const T *a, const T *b, T *c, num_t (*f)(num_t va, num_t vb, int i_))
{
  assert(a && b && c);
  ensure(a->desc == b->desc && a->desc == c->desc);

  idx_t *pi = a->desc->hpoly_To_idx;
  if (a->lo > b->lo) { const T* t; SWAP(a,b,t); }
  ord_t c_hi = min_ord(MAX(a->hi,b->hi),c->mo,c->desc->trunc), c_lo = a->lo;

  num_t va, vb;
  for (int i = pi[c_lo]; i < pi[c_hi+1]; ++i) {
    int curr_strict_mode = mad_tpsa_strict;
    mad_tpsa_strict = 0;
    va = mad_tpsa_geti(a,i);
    vb = mad_tpsa_geti(b,i);
    mad_tpsa_strict = curr_strict_mode;
    c->coef[i] = f(va,vb,i);
  }
  c->lo = c_lo;
  c->hi = c_hi;
  c->nz = mad_bit_trunc(mad_bit_add(a->nz,b->nz),c->hi);

  return c;
}

#undef T
#undef D
#undef TRACE

// --- --- OPERATIONS ---------------------------------------------------------

#include "mad_tpsa_ops.tc"

#include "mad_tpsa_fun.tc"

#include "mad_tpsa_comp.tc"

#include "mad_tpsa_minv.tc"

#include "mad_tpsa_io.tc"
