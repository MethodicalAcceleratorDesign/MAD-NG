/*
 o-----------------------------------------------------------------------------o
 |
 | Monomial module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
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
#include <stdlib.h>
#include <assert.h>

#include "mad_mono.h"

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

ord_t
mad_mono_min (ssz_t n, const ord_t a[n])
{
  assert(a);
  ord_t mo = ~0;
  for (idx_t i=0; i < n; ++i) if (a[i] < mo) mo = a[i];
  return mo;
}

ord_t
mad_mono_max (ssz_t n, const ord_t a[n])
{
  assert(a);
  ord_t mo = 0;
  for (idx_t i=0; i < n; ++i) if (a[i] > mo) mo = a[i];
  return mo;
}

int
mad_mono_ord (ssz_t n, const ord_t a[n])
{
  assert(a);
  int s = 0;
  for (idx_t i=0; i < n; ++i) s += a[i];
  return s;
}

log_t
mad_mono_eq (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] != b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_lt (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] >= b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_le (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] > b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_gt (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] <= b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_ge (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] < b[i]) return FALSE;
  return TRUE;
}

int
mad_mono_cmp (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=0; i < n; ++i) if (a[i] != b[i]) return (int)a[i] - b[i];
  return 0;
}

int
mad_mono_rcmp (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (idx_t i=n-1; i >= 0; --i) if (a[i] != b[i]) return (int)a[i] - b[i];
  return 0;
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
  for (idx_t i=0; i < n; ++i) r[i] = a[i] > b[i] ? a[i] - b[i] : 0;
}

void
mad_mono_cat (ssz_t n, const ord_t a[n],
              ssz_t m, const ord_t b[m], ord_t r[n+m])
{
  mad_mono_copy(n, a, r  );
  mad_mono_copy(m, b, r+n);
}

// -- sorting

static __thread const ord_t *ords;

static int
cmp (const void *a, const void *b)
{
  int i1 = *(const int*)a;
  int i2 = *(const int*)b;

  return (int)ords[i1] - ords[i2];
}

void
mad_mono_sort (ssz_t n, const ord_t a[n], idx_t idxs[n])
{
  assert(a && idxs);
  ords = a;
  for (idx_t i=0; i < n; ++i) idxs[i] = i;
  qsort(idxs, n, sizeof *idxs, cmp);
}

// -- printing

void
mad_mono_print (ssz_t n, const ord_t m[n])
{
  assert(m);

  printf("[ ");
  for (idx_t i=0; i < n; ++i)
    printf("%d ", (int)m[i]);
  printf("]");
}

// --- end --------------------------------------------------------------------o
