/*
 o----------------------------------------------------------------------------o
 |
 | Monomial module implementation
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
#include <assert.h>

#include "mad_mono.h"

// --- locals ----------------------------------------------------------------o

static const ord_t *ords_to_sort;

static int
cmp_ords (const void *a, const void *b)
{
  int i1 = *(const int*)a;
  int i2 = *(const int*)b;

  return ords_to_sort[i1] > ords_to_sort[i2] ? -1 : ords_to_sort[i1] < ords_to_sort[i2];
}

// --- implementation --------------------------------------------------------o

void
mad_mono_fill (int n, ord_t a[n], ord_t v)
{
  assert(a);
  for (int i=0; i < n; ++i) a[i] = v;
}

void
mad_mono_copy (int n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  for (int i=0; i < n; ++i) r[i] = a[i];
}

void
mad_mono_concat (int n, const ord_t a[], int m, const ord_t b[], ord_t r[])
{
  assert(a && b && r);
  int i, j;
  for (i=0; i < n; ++i) r[i  ] = a[i];
  for (j=0; j < m; ++j) r[i+j] = b[i];
}

ord_t
mad_mono_max (int n, const ord_t a[n])
{
  assert(a);
  ord_t mo = a[0];
  for (int i = 1; i < n; ++i)
    if (a[i] > mo) mo = a[i];
  return mo;
}

int
mad_mono_equ (int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (int i=0; i < n; ++i)
    if (a[i] != b[i]) return 0;
  return 1;
}

int
mad_mono_rcmp (int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (int i=n-1; i >= 0; --i)
    if (a[i] != b[i]) return a[i] - b[i];
  return 0;
}

void
mad_mono_sort (int n, const ord_t a[n], int idxs[n])
{
  assert(a && idxs);

  ords_to_sort = a;
  for (int i=0; i < n; ++i) idxs[i] = i;
  qsort(idxs, n, sizeof *idxs, cmp_ords);
}

void
mad_mono_print (int n, const ord_t m[n])
{
  assert(m);

  printf("[ ");
  for (int i=0; i < n; ++i)
    printf("%d ", (int)m[i]);
  printf("]");
}

// --- optimized versions ----------------------------------------------------o

#if   defined(__AVX2__)
#  include "sse/mad_mono_sse.tc"
// NYI -> #  include "avx/mad_mono_avx.tc"
#elif defined(__SSE2__)
#  include "sse/mad_mono_sse.tc"
#else

// --- default versions ------------------------------------------------------o

void
mad_mono_add (int n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (int i=0; i < n; ++i) r[i] = a[i] + b[i];
}

void
mad_mono_sub (int n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (int i=0; i < n; ++i) r[i] = a[i] - b[i];
}

int
mad_mono_leq (int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (int i=0; i < n; ++i)
    if (a[i] > b[i]) return 0;
  return 1;
}

int
mad_mono_ord (int n, const ord_t a[n])
{
  assert(a);
  int s = 0;
  for (int i=0; i < n; ++i) s += a[i];
  return s;
}

#endif // __SSE2__ || __AVX2__
