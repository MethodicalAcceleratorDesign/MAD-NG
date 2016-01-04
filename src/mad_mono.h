#ifndef MAD_MONO_H
#define MAD_MONO_H

// --- types -------------------------------------------------------------------

typedef unsigned char ord_t;

// --- interface ---------------------------------------------------------------

static void  mono_fill    (int n, ord_t a[n], ord_t v);
static ord_t mono_order   (int n, const ord_t a[n]);
static void  mono_copy    (int n, const ord_t a[n],       ord_t r[n]);
static int   mono_equ     (int n, const ord_t a[n], const ord_t b[n]);
static void  mono_add     (int n, const ord_t a[n], const ord_t b[n], ord_t r[n]);
static void  mono_sub     (int n, const ord_t a[n], const ord_t b[n], ord_t r[n]);

// -----------------------------------------------------------------------------
static ord_t mono_max     (int n, const ord_t a[n]);
static int   mono_rcmp    (int n, const ord_t a[n], const ord_t b[n]);
static int   mono_leq     (int n, const ord_t a[n], const ord_t b[n]);
static void  mono_print   (int n, const ord_t a[n]);
static void  mono_sort    (int n, const ord_t a[n], int idxs[n]);
// --- implementation ----------------------------------------------------------

#include <assert.h>

static inline void
mono_fill(int n, ord_t a[n], ord_t v)
{
  assert(a);
  for (int i=0; i < n; ++i) a[i] = v;
}

static inline ord_t
mono_order(int n, const ord_t a[n])
{
  assert(a);
  ord_t s = 0;
  for (int i=0; i < n; ++i) s += a[i];
  return s;
}

static inline ord_t
mono_max(int n, const ord_t a[n])
{
  assert(a);
  ord_t mo = a[0];
  for (int i = 1; i < n; ++i)
    if (a[i] > mo)
      mo = a[i];
  return mo;
}

static inline void
mono_sub(int n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (int i = 0; i < n; ++i) r[i] = a[i] - b[i];
}

static inline void
mono_copy(int n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  for (int i = 0; i < n; ++i) r[i] = a[i];
}

static inline int
mono_equ(int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (int i = 0; i < n; ++i)
    if (a[i] != b[i]) return 0;
  return 1;
}

static inline int
mono_rcmp(int n, const ord_t a[n], const ord_t b[n])
{
  assert(a);
  assert(b);
  for (int i = n - 1; i >= 0; --i)
    if (a[i] != b[i])
      return a[i] - b[i];
  return 0;
}

static inline int
mono_leq(int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  for (int i=0; i < n; ++i)
    if (a[i] > b[i]) return 0;
  return 1;
}

static inline void
mono_add(int n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  for (int i = 0; i < n; ++i) r[i] = a[i] + b[i];
}

static const ord_t *ords_to_sort;

static inline int
cmp_ords(const void *a, const void *b)
{
  int i1 = *(int*)a, i2 = *(int*)b;
  return ords_to_sort[i1] > ords_to_sort[i2] ? -1 : ords_to_sort[i1] < ords_to_sort[i2];
}

static inline void
mono_sort(int n, const ord_t a[n], int idxs[n])
{
  assert(a && idxs);
  ords_to_sort = a;
  for (int i = 0; i < n; ++i)
    idxs[i] = i;
  qsort(idxs,n,sizeof(*idxs),cmp_ords);
}

#include <stdio.h>

static inline void
mono_print(int n, const ord_t m[n])
{
  assert(m);
  printf("[ ");
  for (int i=0; i < n; ++i)
    printf("%d ", (int)m[i]);
  printf("]");
}

// --- SSE2 implementation -----------------------------------------------------

// Comment the following include to disable SSE/AVX optimization
#include "mono_sse.h"
#include "mono_avx.h"

// -----------------------------------------------------------------------------
#endif
