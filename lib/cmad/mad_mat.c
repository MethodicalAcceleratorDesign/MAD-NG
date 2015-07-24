#include "mad_mat.h"
#include <complex.h>
#include <assert.h>

typedef double           num_t;
typedef double _Complex cnum_t;

#define CHKXY  assert( x && y )
#define CHKXR  assert( x && y && r && x != r )
#define CHKYR  assert( x && y && r && y != r )
#define CHKXYR assert( x && y && r && x != r && y != r )

// [m x n] = [m x p] * [p x n]
/* Naive implementation (not vectorized)
#define MMUL \
  for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
*/

// portable vectorized matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
// get xN speed-up factor compared to dgemm from openblas and lapack...

// [m x n] = [m x p] * [p x n]
#define MMUL { \
  if (n >= 8) { \
    for (size_t i=0; i < m; i++) { \
      for (size_t j=0; j < n-7; j+=8) { \
        r[i*n+j  ] = r[i*n+j+1] = \
        r[i*n+j+2] = r[i*n+j+3] = \
        r[i*n+j+4] = r[i*n+j+5] = \
        r[i*n+j+6] = r[i*n+j+7] = 0; \
        for (size_t k=0; k < p; k++) { \
          r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
          r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
          r[i*n+j+2] += x[i*p+k] * y[k*n+j+2]; \
          r[i*n+j+3] += x[i*p+k] * y[k*n+j+3]; \
          r[i*n+j+4] += x[i*p+k] * y[k*n+j+4]; \
          r[i*n+j+5] += x[i*p+k] * y[k*n+j+5]; \
          r[i*n+j+6] += x[i*p+k] * y[k*n+j+6]; \
          r[i*n+j+7] += x[i*p+k] * y[k*n+j+7]; \
        } \
      } \
    } \
  } \
  if (n & 4) { \
    size_t j = n - n%8; \
    for (size_t i=0; i < m; i++) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = 0; \
      for (size_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
        r[i*n+j+2] += x[i*p+k] * y[k*n+j+2]; \
        r[i*n+j+3] += x[i*p+k] * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    size_t j = n - n%4; \
    for (size_t i=0; i < m; i++) { \
      r[i*n+j] = r[i*n+j+1] = 0; \
      for (size_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    size_t j = n - 1; \
    for (size_t i=0; i < m; i++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
  } \
} \

// n==1: [m x 1] = [m x p] * [p x 1]
#define MULV \
  for (size_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (size_t k=0; k < p; k++) \
      r[i] += x[i*p+k] * y[k]; \
  } \

// m==1: [1 x n] = [1 x p] * [p x n]
#define VMUL \
  for (size_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (size_t k=0; k < p; k++) \
      r[j] += x[k] * y[k*n+j]; \
  } \

// p==1: [m x n] = [m x 1] * [1 x n]
#define KMUL \
  for (size_t i=0; i < m; i++) { \
  for (size_t j=0; j < n; j++) \
    r[i*n+j] = x[i] * y[j]; \
  } \

// m==1, n==1: [1 x 1] = [1 x p] * [p x 1]
#define VMULV \
  { r[0] = 0; \
    for (size_t k=0; k < p; k++) \
      r[0] += x[k] * y[k]; \
  } \

// -----

#define MUL() { \
  if (m == 1 && n == 1) VMULV \
  else      if (m == 1) VMUL \
  else      if (n == 1) MULV \
  else      if (p == 1) KMUL \
  else                  MMUL \
} \

#define DOT(C) { \
  *r = 0; \
  if (m == 1 && n == 1) { \
    for (size_t i=0; i < p; i++) \
      *r += C(x[i]) * y[i]; \
  } else { \
    size_t mn = m < n ? m : n; \
    for (size_t i=0; i < mn; i++) \
      for (size_t k=0; k < p; k++) \
        *r += C(x[k*m+i]) * y[k*n+i]; \
  } \
} \

// transpose, in place only for square matrix
// inplace for non-square matrix could use FFTW...
// see http://www.fftw.org/faq/section3.html#transpose
#define TRANS(C) { \
  assert(x && r); \
  if (m == 1 || n == 1) { \
    size_t mn = m*n; \
    for (size_t i=0; i < mn; i++) \
      r[i] = C(x[i]); \
  } else if (x == r) { \
    assert(m == n); \
    for (size_t i=0; i < m; i++) \
    for (size_t j=i; j < n; j++) \
      r[j*n+i] = C(r[i*n+j]); \
  } else { \
    for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) \
      r[j*m+i] = C(x[i*n+j]); \
  } \
} \


void mad_mat_trans (const num_t *x, num_t *r, size_t m, size_t n)
{ TRANS(); }

num_t mad_mat_dot (const num_t *x, const num_t *y, size_t m, size_t n, size_t p)
{ CHKXY; num_t r_, *r=&r_; DOT(); return *r; }

void mad_mat_dotm (const num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(); }

void mad_mat_mul (const num_t *x, const num_t *y, num_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; MUL(); }

void mad_mat_mulm (const num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKYR; MUL(); }

void mad_cmat_trans (const cnum_t *x, cnum_t *r, size_t m, size_t n)
{ TRANS(); }

void mad_cmat_ctrans (const cnum_t *x, cnum_t *r, size_t m, size_t n)
{ TRANS(conj); }

void mad_cmat_dot (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(conj); }

void mad_cmat_dotm (const cnum_t *x, const num_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(conj); }

void mad_cmat_mul (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; MUL(); }

void mad_cmat_mulm (const cnum_t *x, const num_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXR; MUL(); }
