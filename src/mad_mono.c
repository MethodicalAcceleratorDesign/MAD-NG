/*
 o-----------------------------------------------------------------------------o
 |
 | Monomial module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: C. Tomoiaga
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mad_mono.h"

// --- locals -----------------------------------------------------------------o

static const ord_t *ords_to_sort;

static int
cmp_ords (const void *a, const void *b)
{
  int i1 = *(const int*)a;
  int i2 = *(const int*)b;

  return ords_to_sort[i1] > ords_to_sort[i2] ? -1 :
         ords_to_sort[i1] < ords_to_sort[i2];
}

// --- implementation ---------------------------------------------------------o

ssz_t
mad_mono_str (ssz_t n, ord_t a[n], str_t s)
{
  assert(a && s);
  ssz_t i = 0;
  for (; i < n && s[i]; ++i) a[i] = s[i] - (s[i] > '9' ? 'A'-10 : '0');
  return i;
}

void
mad_mono_fill (ssz_t n, ord_t a[n], ord_t v)
{
  assert(a);
  for (idx_t i=0; i < n; ++i) a[i] = v;
}

void
mad_mono_copy (ssz_t n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  if (a != r)
    for (idx_t i=0; i < n; ++i) r[i] = a[i];
}

void
mad_mono_rcopy (ssz_t n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  ord_t t;
  if (a != r)
    for (idx_t i=0; i < n; ++i) r[i] = a[n-1-i];
  else
    for (idx_t i=0; i < n/2; ++i)
      t = r[i], r[i] = r[n-1-i], r[n-1-i] = t;
}

void
mad_mono_add (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (idx_t i=0; i < n; ++i) r[i] = a[i] + b[i];
}

void
mad_mono_sub (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (idx_t i=0; i < n; ++i) r[i] = a[i] - b[i];
}

int
mad_mono_ord (ssz_t n, const ord_t a[n])
{
  assert(a);
  int s = 0;
  for (idx_t i=0; i < n; ++i) s += a[i];
  return s;
}

void
mad_mono_concat (ssz_t n, const ord_t a[n],
                 ssz_t m, const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  idx_t i, j;
  for (i=0; i < n; ++i) r[i  ] = a[i];
  for (j=0; j < m; ++j) r[i+j] = b[j];
}

void
mad_mono_sort (ssz_t n, const ord_t a[n], idx_t idxs[n])
{
  assert(a && idxs);
  ords_to_sort = a;
  for (idx_t i=0; i < n; ++i) idxs[i] = i;
  qsort(idxs, n, sizeof *idxs, cmp_ords);
}

void
mad_mono_print (ssz_t n, const ord_t m[n])
{
  assert(m);

  printf("[ ");
  for (idx_t i=0; i < n; ++i)
    printf("%d ", (int)m[i]);
  printf("]");
}

// --- default/reference versions ---------------------------------------------o

#undef __SSE2__ // disable optimized versions

#ifdef __SSE2__
#define FUN(name)   MKNAME(name,_ref)
#else
#define FUN(name)   name
#endif

int
FUN(mad_mono_eq) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i)
    if (a[i] != b[i]) return 0;
  return 1;
}

int
FUN(mad_mono_lt) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i)
    if (a[i] >= b[i]) return 0;
  return 1;
}

int
FUN(mad_mono_gt) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i)
    if (a[i] <= b[i]) return 0;
  return 1;
}

int
FUN(mad_mono_le) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  return !FUN(mad_mono_gt)(n,a,b);
}

int
FUN(mad_mono_ge) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  return !FUN(mad_mono_lt)(n,a,b);
}

int
FUN(mad_mono_rcmp) (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=n-1; i >= 0; --i)
    if (a[i] - b[i]) return (int)a[i] - b[i];
  return 0;
}

ord_t
FUN(mad_mono_max) (ssz_t n, const ord_t a[n])
{
  assert(a);
  ord_t mo = 0;
  for (idx_t i=0; i < n; ++i)
    if (a[i] > mo) mo = a[i];
  return mo;
}

ord_t
FUN(mad_mono_min) (ssz_t n, const ord_t a[n])
{
  assert(a);
  ord_t mo = ~0;
  for (idx_t i=0; i < n; ++i)
    if (a[i] < mo) mo = a[i];
  return mo;
}

// --- optimized versions -----------------------------------------------------o

// to check availability with gcc:
// echo | gcc -dM -E -march=native - | grep "SSE\|AVX"
// echo | gcc -dM -E -msse2 - | grep "SSE\|AVX"
// echo | gcc -dM -E -mavx2 - | grep "SSE\|AVX"

#if 0 // disable optimized versions

#if defined(__AVX512F__) && defined(__AVX512BW__)
// #warning "AVX512 selected"
#include "sse/mad_mono_avx512.tc" // never tested
#elif defined(__AVX2__)
// #warning "AVX2 selected"
#include "sse/mad_mono_avx2.tc"
#elif defined(__SSE2__)
// #warning "SSE2 selected"
#include "sse/mad_mono_sse2.tc"
#endif // __SSE2__ || __AVX2__ || __AVX512F__

#endif

// --- end --------------------------------------------------------------------o
