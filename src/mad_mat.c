/*
 o-----------------------------------------------------------------------------o
 |
 | Matrix module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"

// --- implementation ---------------------------------------------------------o

#define CHKR     assert( r )
#define CHKX     assert( x )
#define CHKXY    assert( x && y )
#define CHKXR    assert( x && r )
#define CHKYR    assert( y && r )
#define CHKXYR   assert( x && y && r )
#define CHKXRX   assert( x && r && x != r)

#define CNUM(a) cnum_t a = (* (cnum_t*) & (num_t[2]) { MKNAME(a,_re), MKNAME(a,_im) })

// -----

#if 0
// [m x n] = [m x p] * [p x n]
// naive implementation (not vectorized)
#define MMUL { /* mat * mat */ \
  for (ssz_t i=0; i < m; i++) \
  for (ssz_t j=0; j < n; j++) { \
    r[i*n+j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i*n+j] += x[i*p+k] * y[k*n+j]; \
  } \
}
#else
// [m x n] = [m x p] * [p x n]
// portable vectorized general matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define MMUL() { /* mat * mat */ \
  assert(m>0 && n>0 && p>0); \
  if (n & ~7) { \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n-7; j+=8) { \
        r[i*n+j  ] = r[i*n+j+1] = \
        r[i*n+j+2] = r[i*n+j+3] = \
        r[i*n+j+4] = r[i*n+j+5] = \
        r[i*n+j+6] = r[i*n+j+7] = 0; \
        for (ssz_t k=0; k < p; k++) { \
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
    ssz_t j = n - (n & 7); \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = 0; \
      for (ssz_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
        r[i*n+j+2] += x[i*p+k] * y[k*n+j+2]; \
        r[i*n+j+3] += x[i*p+k] * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    ssz_t j = n - (n & 3); \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j] = r[i*n+j+1] = 0; \
      for (ssz_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    ssz_t j = n - 1; \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j] = 0; \
      for (ssz_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
  } \
}
#endif

// n=1: [m x 1] = [m x p] * [p x 1]
#define MULV() /* mat * vec */ \
  for (ssz_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i] += x[i*p+k] * y[k]; \
  }

// m=1: [1 x n] = [1 x p] * [p x n]
#define VMUL() /* vec * mat */ \
  for (ssz_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[j] += x[k] * y[k*n+j]; \
  }

// m=1, n=1: [1 x 1] = [1 x p] * [p x 1]
#define IMUL() /* vec * vec */ \
  { r[0] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[0] += x[k] * y[k]; \
  }

// [m x n] = [m x p] * [p x n]
#define MUL() \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: MMUL(); break; \
    case 1: MULV(); break; \
    case 2: VMUL(); break; \
    case 3: IMUL(); break; \
  }

// -----

#if 0
// [m x n] = [p x m]' * [p x n]
// naive implementation (not vectorized)
#define TMMUL(C) { /* mat' * mat */ \
  for (ssz_t i=0; i < m; i++) \
  for (ssz_t j=0; j < n; j++) { \
    r[i*n+j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i*n+j] += C(x[k*m+i]) * y[k*n+j]; \
  } \
}
#else
// [m x n] = [p x m]' * [p x n]
// portable vectorized general transpose matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define TMMUL(C) { /* mat' * mat */ \
  assert(m>0 && n>0 && p>0); \
  if (n & ~7) { \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n-7; j+=8) { \
        r[i*n+j  ] = r[i*n+j+1] = \
        r[i*n+j+2] = r[i*n+j+3] = \
        r[i*n+j+4] = r[i*n+j+5] = \
        r[i*n+j+6] = r[i*n+j+7] = 0; \
        for (ssz_t k=0; k < p; k++) { \
          r[i*n+j  ] += C(x[k*m+i]) * y[k*n+j  ]; \
          r[i*n+j+1] += C(x[k*m+i]) * y[k*n+j+1]; \
          r[i*n+j+2] += C(x[k*m+i]) * y[k*n+j+2]; \
          r[i*n+j+3] += C(x[k*m+i]) * y[k*n+j+3]; \
          r[i*n+j+4] += C(x[k*m+i]) * y[k*n+j+4]; \
          r[i*n+j+5] += C(x[k*m+i]) * y[k*n+j+5]; \
          r[i*n+j+6] += C(x[k*m+i]) * y[k*n+j+6]; \
          r[i*n+j+7] += C(x[k*m+i]) * y[k*n+j+7]; \
        } \
      } \
    } \
  } \
  if (n & 4) { \
    ssz_t j = n - (n & 7); \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = 0; \
      for (ssz_t k=0; k < p; k++) { \
        r[i*n+j  ] += C(x[k*m+i]) * y[k*n+j  ]; \
        r[i*n+j+1] += C(x[k*m+i]) * y[k*n+j+1]; \
        r[i*n+j+2] += C(x[k*m+i]) * y[k*n+j+2]; \
        r[i*n+j+3] += C(x[k*m+i]) * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    ssz_t j = n - (n & 3); \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j] = r[i*n+j+1] = 0; \
      for (ssz_t k=0; k < p; k++) { \
        r[i*n+j  ] += C(x[k*m+i]) * y[k*n+j  ]; \
        r[i*n+j+1] += C(x[k*m+i]) * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    ssz_t j = n - 1; \
    for (ssz_t i=0; i < m; i++) { \
      r[i*n+j] = 0; \
      for (ssz_t k=0; k < p; k++) \
        r[i*n+j] += C(x[k*m+i]) * y[k*n+j]; \
    } \
  } \
}
#endif

// n=1: [m x 1] = [p x m]' * [p x 1]
#define TMULV(C) /* mat' * vec */ \
  for (ssz_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i] += C(x[k*m+i]) * y[k]; \
  }

// m=1: [1 x n] = [p x 1]' * [p x n]
#define TVMUL(C) /* vec' * mat */ \
  for (ssz_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[j] += C(x[k]) * y[k*n+j]; \
  }

// m=1, n=1: [1 x 1] = [p x 1]' * [p x 1]
#define TIMUL(C) /* vec' * vec */ \
  { r[0] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[0] += C(x[k]) * y[k]; \
  }

// [m x n] = [p x m]' * [p x n]
#define TMUL(C) \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: TMMUL(C); break; \
    case 1: TMULV(C); break; \
    case 2: TVMUL(C); break; \
    case 3: TIMUL(C); break; \
  }

// -----

#if 0
// [m x n] = [m x p] * [n x p]'
// naive implementation (not vectorized)
#define MMULT(C) { /* mat * mat' */ \
  for (ssz_t i=0; i < m; i++) \
  for (ssz_t j=0; j < n; j++) { \
    r[i*n+j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i*n+j] += x[i*p+k] * C(y[j*p+k]); \
  } \
}
#else
// [m x n] = [m x p] * [n x p]'
// portable vectorized general transpose matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define MMULT(C) { /* mat * mat' */ \
  assert(m>0 && n>0 && p>0); \
  for (ssz_t i=0; i < m*n; i++) r[i] = 0; \
  if (n & ~7) { \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n; j++) { \
        for (ssz_t k=0; k < p-7; k+=8) { \
          r[i*n+j] += x[i*p+k  ] * C(y[j*p+k  ]); \
          r[i*n+j] += x[i*p+k+1] * C(y[j*p+k+1]); \
          r[i*n+j] += x[i*p+k+2] * C(y[j*p+k+2]); \
          r[i*n+j] += x[i*p+k+3] * C(y[j*p+k+3]); \
          r[i*n+j] += x[i*p+k+4] * C(y[j*p+k+4]); \
          r[i*n+j] += x[i*p+k+5] * C(y[j*p+k+5]); \
          r[i*n+j] += x[i*p+k+6] * C(y[j*p+k+6]); \
          r[i*n+j] += x[i*p+k+7] * C(y[j*p+k+7]); \
        } \
      } \
    } \
  } \
  if (p & 4) { \
    ssz_t k = p - (p & 7); \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n; j++) { \
        r[i*n+j] += x[i*p+k  ] * C(y[j*p+k  ]); \
        r[i*n+j] += x[i*p+k+1] * C(y[j*p+k+1]); \
        r[i*n+j] += x[i*p+k+2] * C(y[j*p+k+2]); \
        r[i*n+j] += x[i*p+k+3] * C(y[j*p+k+3]); \
      } \
    } \
  } \
  if (p & 2) { \
    ssz_t k = p - (p & 3); \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n; j++) { \
        r[i*n+j] += x[i*p+k  ] * C(y[j*p+k  ]); \
        r[i*n+j] += x[i*p+k+1] * C(y[j*p+k+1]); \
      } \
    } \
  } \
  if (p & 1) { \
    ssz_t k = p - 1; \
    for (ssz_t i=0; i < m; i++) { \
      for (ssz_t j=0; j < n; j++) \
        r[i*n+j] += x[i*p+k] * C(y[j*p+k]); \
    } \
  } \
}
#endif

// n=1: [m x 1] = [m x p] * [1 x p]'
#define MULVT(C) /* mat * vec' */ \
  for (ssz_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[i] += x[i*p+k] * C(y[k]); \
  }

// m=1: [1 x n] = [1 x p]' * [n x p]'
#define VMULT(C) /* vec * mat' */ \
  for (ssz_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[j] += x[k] * C(y[j*p+k]); \
  }

// m=1, n=1: [1 x 1] = [1 x p] * [1 x p]'
#define IMULT(C) /* vec * vec' */ \
  { r[0] = 0; \
    for (ssz_t k=0; k < p; k++) \
      r[0] += x[k] * C(y[k]); \
  }

// [m x n] = [m x p] * [n x p]'
#define MULT(C) \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: MMULT(C); break; \
    case 1: MULVT(C); break; \
    case 2: VMULT(C); break; \
    case 3: IMULT(C); break; \
  }

// -----

// r = ([m x n] <*> [m x n])
#define DOT(C) { \
  if (n == 1) { \
    r[0] = 0; \
    for (ssz_t i=0; i < m; i++) \
      r[0] += C(x[i]) * y[i]; \
  } else { \
    for (ssz_t j=0; j < n; j++) { \
      r[j] = 0; \
      for (ssz_t i=0; i < m; i++) \
        r[j] += C(x[i*n+j]) * y[i*n+j]; \
    } \
  } \
}

// [m x n] transpose
#define TRANS(C,T) { \
  if (m == 1 || n == 1) { \
    ssz_t mn = m*n; \
    for (ssz_t i=0; i < mn; i++) \
      r[i] = C(x[i]); \
  } else if ((const void*)x != (const void*)r) { \
    for (ssz_t i=0; i < m; i++) \
    for (ssz_t j=0; j < n; j++) \
      r[j*m+i] = C(x[i*n+j]); \
  } else if (m == n) { \
    for (ssz_t i=0; i < m; i++) \
    for (ssz_t j=i; j < n; j++) { \
      T t = C(r[j*m+i]); \
      r[j*m+i] = C(r[i*n+j]); \
      r[i*n+j] = t; \
    } \
  } else { \
    mad_alloc_tmp(T, t, m*n); \
    for (ssz_t i=0; i < m; i++) \
    for (ssz_t j=0; j < n; j++) \
      t[j*m+i] = C(x[i*n+j]); \
    memcpy(r, t, m*n*sizeof(T)); \
    mad_free_tmp(t); \
  } \
}

#define CPY(OP) { \
  for (ssz_t i=0; i<m; i++) \
  for (ssz_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x[i*ldx+j]; \
}

#define SET(OP) { \
  for (ssz_t i=0; i<m; i++) \
  for (ssz_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x; \
}

#define DIAG(OP) { \
  for (ssz_t i=0; i<MIN(m,n); i++) \
    r[i*ldr+i] OP##= x; \
}

// --- mat

void mad_mat_reshape(struct matrix *x, ssz_t m, ssz_t n)
{ assert(x); x->nr = m; x->nc = n; }

void mad_mat_ident(num_t r[], ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; num_t x = 0; SET(); x = 1; DIAG(); }

void mad_mat_fill(num_t x, num_t r[], ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; SET(); }

void mad_mat_copy(const num_t x[], num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXRX; CPY(); }

void mad_mat_copym(const num_t x[], cnum_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXR; CPY(); }

void mad_mat_trans (const num_t x[], num_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(,num_t); }

void mad_mat_dot (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n)
{ CHKXYR; DOT(); }

void mad_mat_dotm (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n)
{ CHKXYR; DOT(); }

void mad_mat_mul (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MUL(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  MUL();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_mulm (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { MUL(); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_tmul (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { TMUL(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  TMUL();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_tmulm (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { TMUL(); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  TMUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_mult (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MULT(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  MULT();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_multm (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { MULT(conj); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MULT(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_center (const num_t x[], num_t r[], ssz_t m, ssz_t n, int d)
{ CHKXR;
  assert(d == 1 || d == 2); // 1=row, 2=col
  if (d == 1)
    for (ssz_t i=0; i < m; i++) {
      num_t mu = 0;
      for (ssz_t j=0; j < n; j++) mu += x[i*n+j];
      mu /= n;
      for (ssz_t j=0; j < n; j++) r[i*n+j] = x[i*n+j] - mu;
    }
  else
    for (ssz_t j=0; j < n; j++) {
      num_t mu = 0;
      for (ssz_t i=0; i < m; i++) mu += x[i*n+j];
      mu /= m;
      for (ssz_t i=0; i < n; i++) r[i*n+j] = x[i*n+j] - mu;
    }
}

void
mad_mat_roll (num_t x[], ssz_t m, ssz_t n, int mroll, int nroll)
{ CHKX; mroll %= m; nroll %= n;
  ssz_t nm = n*m, msz = n*abs(mroll), nsz = abs(nroll);
  ssz_t sz = msz > nsz ? msz : nsz;
  mad_alloc_tmp(num_t, a, sz);
  if (mroll > 0) {
    mad_vec_copy(x+nm-msz, a    ,    msz); // end of x to a
    mad_vec_copy(x       , x+msz, nm-msz); // shift x down
    mad_vec_copy(a       , x    ,    msz); // a to beginning of x
  } else
  if (mroll < 0) {
    mad_vec_copy(x    , a       ,    msz); // beginning of x to a
    mad_vec_copy(x+msz, x       , nm-msz); // shift x up
    mad_vec_copy(a    , x+nm-msz,    msz); // a to end of x
  }
  if (nroll > 0) {
    for (ssz_t i=0; i < nm; i += n) {
      mad_vec_copy(x+i+n-nsz, a      ,   nsz); // end of x to a
      mad_vec_copy(x+i      , x+i+nsz, n-nsz); // shift x right
      mad_vec_copy(a        , x+i    ,   nsz); // a to beginning of x
    }
  } else
  if (nroll < 0) {
    for (ssz_t i=0; i < nm; i += n) {
      mad_vec_copy(x+i    , a        ,   nsz); // beginning of x to a
      mad_vec_copy(x+i+nsz, x+i      , n-nsz); // shift x left
      mad_vec_copy(a      , x+i+n-nsz,   nsz); // a to end of x
    }
  }
  mad_free_tmp(a);
}

// -- cmat

void mad_cmat_reshape(struct cmatrix *x, ssz_t m, ssz_t n)
{ assert(x); x->nr = m; x->nc = n; }

void mad_cmat_ident(cnum_t r[], ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; cnum_t x = 0; SET(); x = 1; DIAG(); }

void mad_cmat_fill(cnum_t x, cnum_t r[], ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; SET(); }

void mad_cmat_fill_r(num_t x_re, num_t x_im, cnum_t r[], ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; CNUM(x); SET(); }

void mad_cmat_roll (cnum_t x[], ssz_t m, ssz_t n, int mroll, int nroll)
{ mad_mat_roll((num_t*)x, m, 2*n, mroll, 2*nroll); }

void mad_cmat_copy(const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXRX; CPY(); }

void mad_cmat_trans (const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(,cnum_t); }

void mad_cmat_ctrans (const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(conj,cnum_t); }

void mad_cmat_dot (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n)
{ CHKXYR; DOT(conj); }

void mad_cmat_dotm (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t m, ssz_t n)
{ CHKXYR; DOT(conj); }

void mad_cmat_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MUL(); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_mulm (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { MUL(); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_tmul (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { TMUL(conj); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  TMUL(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_tmulm (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { TMUL(conj); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  TMUL(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_mult (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MULT(conj); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MULT(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_multm (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { MULT(); return; }
  mad_alloc_tmp(cnum_t, r_, m*n);
  cnum_t *t = r; r = r_;
  MULT();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_center (const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n, int d)
{ CHKXR;
  assert(d == 1 || d == 2); // 1=row, 2=col
  if (d == 1)
    for (ssz_t i=0; i < m; i++) {
      cnum_t mu = 0;
      for (ssz_t j=0; j < n; j++) mu += x[i*n+j];
      mu /= n;
      for (ssz_t j=0; j < n; j++) r[i*n+j] = x[i*n+j] - mu;
    }
  else
    for (ssz_t j=0; j < n; j++) {
      cnum_t mu = 0;
      for (ssz_t i=0; i < m; i++) mu += x[i*n+j];
      mu /= m;
      for (ssz_t i=0; i < n; i++) r[i*n+j] = x[i*n+j] - mu;
    }
}

// -- Symplectic matrices -----------------------------------------------------o

// M[2n x 2n] accessed as n blocks of [a b ; c d]

#define a_(x,i,j) x[ i   *n+j  ]
#define b_(x,i,j) x[ i   *n+j+1]
#define c_(x,i,j) x[(i+1)*n+j  ]
#define d_(x,i,j) x[(i+1)*n+j+1]

// -- Symplecticity error, compute M' J M - J ---------------------------------o

num_t mad_mat_symperr (const num_t x[], num_t r[], ssz_t n)
{ CHKX; assert(x != r && !(n & 1));
  num_t s=0, s0, s1, s2, s3;
  for (ssz_t i = 0; i < n-1; i += 2) {
    // i == j
    s1 = -1, s2 = 1;
    for (ssz_t k = 0; k < n-1; k += 2) {
      s1 += a_(x,k,i) * d_(x,k,i) - b_(x,k,i) * c_(x,k,i);
      s2 += b_(x,k,i) * c_(x,k,i) - a_(x,k,i) * d_(x,k,i);
    }
    s += s1*s1 + s2*s2;
    if (r) b_(r,i,i) = s1, c_(r,i,i) = s2, a_(r,i,i) = d_(r,i,i) = 0;
    // i < j
    for (ssz_t j = i+2; j < n-1; j += 2) {
      s0 = s1 = s2 = s3 = 0;
      for (ssz_t k = 0; k < n-1; k += 2) {
        s0 += a_(x,k,i) * c_(x,k,j) - a_(x,k,j) * c_(x,k,i);
        s1 += a_(x,k,i) * d_(x,k,j) - b_(x,k,j) * c_(x,k,i);
        s2 += b_(x,k,i) * c_(x,k,j) - a_(x,k,j) * d_(x,k,i);
        s3 += b_(x,k,i) * d_(x,k,j) - b_(x,k,j) * d_(x,k,i);
      }
      s += 2*(s0*s0 + s1*s1 + s2*s2 + s3*s3);
      if (r) {
        a_(r,i,j) =  s0, b_(r,i,j) =  s1, c_(r,i,j) =  s2, d_(r,i,j) =  s3;
        a_(r,j,i) = -s0, b_(r,j,i) = -s2, c_(r,j,i) = -s1, d_(r,j,i) = -s3;
      }
    }
  }
  return sqrt(s);
}

num_t mad_cmat_symperr (const cnum_t x[], cnum_t r[], ssz_t n)
{ CHKX; assert(x != r && !(n & 1));
  cnum_t s=0, s0, s1, s2, s3;
  for (ssz_t i = 0; i < n-1; i += 2) {
    // i == j
    s1 = -1, s2 = 1;
    for (ssz_t k = 0; k < n-1; k += 2) {
      s1 += conj(a_(x,k,i)) * d_(x,k,i) - b_(x,k,i) * conj(c_(x,k,i));
      s2 += conj(b_(x,k,i)) * c_(x,k,i) - a_(x,k,i) * conj(d_(x,k,i));
    }
    s += s1*s1 + s2*s2;
    if (r) b_(r,i,i) = s1, c_(r,i,i) = s2, a_(r,i,i) = d_(r,i,i) = 0;
    // i < j
    for (ssz_t j = i+2; j < n-1; j += 2) {
      s0 = s1 = s2 = s3 = 0;
      for (ssz_t k = 0; k < n-1; k += 2) {
        s0 += conj(a_(x,k,i)) * c_(x,k,j) - a_(x,k,j) * conj(c_(x,k,i));
        s1 += conj(a_(x,k,i)) * d_(x,k,j) - b_(x,k,j) * conj(c_(x,k,i));
        s2 += conj(b_(x,k,i)) * c_(x,k,j) - a_(x,k,j) * conj(d_(x,k,i));
        s3 += conj(b_(x,k,i)) * d_(x,k,j) - b_(x,k,j) * conj(d_(x,k,i));
      }
      s += 2*(s0*s0 + s1*s1 + s2*s2 + s3*s3);
      if (r) {
        a_(r,i,j) =  s0, b_(r,i,j) =  s1, c_(r,i,j) =  s2, d_(r,i,j) =  s3;
        a_(r,j,i) = -s0, b_(r,j,i) = -s2, c_(r,j,i) = -s1, d_(r,j,i) = -s3;
      }
    }
  }
  return sqrt(cabs(s));
}

// -- Symplectic conjugate, compute \bar{M} = -J M' J -------------------------o

void mad_mat_sympconj (const num_t x[], num_t r[], ssz_t n)
{ CHKXR; assert(!(n & 1));
  num_t t;
  for (ssz_t i = 0; i < n-1; i += 2) {     // 2x2 blocks on diagonal
    t = a_(x,i,i),  a_(r,i,i) =  d_(x,i,i),  d_(r,i,i) = t;
    b_(r,i,i) = -b_(x,i,i),  c_(r,i,i) = -c_(x,i,i);

    for (ssz_t j = i+2; j < n-1; j += 2) { // 2x2 blocks off diagonal
      t = a_(x,i,j),  a_(r,i,j) =  d_(x,j,i),  d_(r,j,i) =  t;
      t = b_(x,i,j),  b_(r,i,j) = -b_(x,j,i),  b_(r,j,i) = -t;
      t = c_(x,i,j),  c_(r,i,j) = -c_(x,j,i),  c_(r,j,i) = -t;
      t = d_(x,i,j),  d_(r,i,j) =  a_(x,j,i),  a_(r,j,i) =  t;
    }
  }
}

void mad_cmat_sympconj (const cnum_t x[], cnum_t r[], ssz_t n)
{ CHKXR; assert(!(n & 1));
  cnum_t t;
  for (ssz_t i = 0; i < n-1; i += 2) {     // 2x2 blocks on diagonal
    t = a_(x,i,i),  a_(r,i,i) =  conj(d_(x,i,i)),  d_(r,i,i) = conj(t);
    b_(r,i,i) = -conj(b_(x,i,i)),  c_(r,i,i) = -conj(c_(x,i,i));

    for (ssz_t j = i+2; j < n-1; j += 2) {   // 2x2 blocks off diagonal
      t = a_(x,i,j),  a_(r,i,j) =  conj(d_(x,j,i)),  d_(r,j,i) =  conj(t);
      t = b_(x,i,j),  b_(r,i,j) = -conj(b_(x,j,i)),  b_(r,j,i) = -conj(t);
      t = c_(x,i,j),  c_(r,i,j) = -conj(c_(x,j,i)),  c_(r,j,i) = -conj(t);
      t = d_(x,i,j),  d_(r,i,j) =  conj(a_(x,j,i)),  a_(r,j,i) =  conj(t);
    }
  }
}

#undef a_
#undef b_
#undef c_
#undef d_

// -- lapack ------------------------------------------------------------------o

/*
LAPACK is the default method for computing LU decomposition. When matrix is
square and nonsinguler the routines dgetrf and zgetrf.

LAPACK is the default method for solving dense numerical matrices. When the
matrix is square and nonsingular the routines dgesv and zgesv are used
otherwise routines dgelsy and zgelsy are used.

LAPACK is the default method for computing the entire set of singluar values
and singular vectors. For generalized SVD the routines dgesdd and zgesdd are
used.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. For simple eigenvalues the routines dgeev and zgeev are used. For
generalized eigenvalues the routines dggev and zggev are used.
*/

// -----
// Decompose A = LU with A[m x n] (generalized)
// -----
void dgetrf_ (const int *m, const int *n,  num_t A[], const int *lda,
              int *IPIV, int *info);
void zgetrf_ (const int *m, const int *n, cnum_t A[], const int *lda,
              int *IPIV, int *info);

// -----
// Solve A * X = B with A[n x n], B[n x nrhs] and X[n x nrhs]: search min | b - Ax | using LU
// -----
void dgesv_ (const int *n, const int *nrhs, num_t  A[], const int *lda,
                                 int *IPIV, num_t  B[], const int *ldb, int *info);
void zgesv_ (const int *n, const int *nrhs, cnum_t A[], const int *lda,
                                 int *IPIV, cnum_t B[], const int *ldb, int *info);

// -----
// Solve A * X = B with A[m x n], B[m x nrhs] and X[m x nrhs]: search min | b - Ax | using QR
// -----
void dgelsy_ (const int *m, const int *n, const int *nrhs,
              num_t A[], const int *lda, num_t B[], const int *ldb,
              int jpvt[], const num_t *rcond, int *rank,
              num_t work[], const int lwork[], int *info);
void zgelsy_ (const int *m, const int *n, const int *nrhs,
              cnum_t A[], const int *lda, cnum_t B[], const int *ldb,
              int jpvt[], const num_t *rcond, int *rank,
              cnum_t work[], const int lwork[], num_t rwork[], int *info);

// -----
// SVD A[m x n]
// -----
void dgesdd_ (str_t jobz, const int *m, const int *n, num_t A[], const int *lda,
              num_t S[], num_t U[], const int *ldu, num_t VT[], const int *ldvt,
              num_t work[], int *lwork, int iwork[], int *info);
void zgesdd_ (str_t jobz, const int *m, const int *n, cnum_t A[], const int *lda,
              num_t S[], cnum_t U[], const int *ldu, cnum_t VT[], const int *ldvt,
              cnum_t work[], int *lwork, num_t rwork[], int iwork[], int *info);

// -----
// Eigen values/vectors A[n x n]
// -----
void dgeev_ (str_t jobvl, str_t jobvr, const int *n, num_t A[], const int *lda,
             num_t WR[], num_t WI[],
             num_t VL[], const int *ldvl, num_t VR[], const int *ldvr,
             num_t work[], int *lwork, int *info);
void zgeev_ (str_t jobvl, str_t jobvr, const int *n, cnum_t A[], const int *lda,
             cnum_t W[], cnum_t VL[], const int *ldvl, cnum_t VR[], const int *ldvr,
             cnum_t work[], int *lwork, num_t rwork[], int *info);

// -- determinant -------------------------------------------------------------o

int
mad_mat_det (const num_t x[], num_t *r, ssz_t n)
{
  CHKX;
  const int nn=n;
  int info=0, ipiv[n];
  mad_alloc_tmp(num_t, a, n*n);
  mad_vec_copy(x, a, n*n);
  dgetrf_(&nn, &nn, a, &nn, ipiv, &info);

  if (info < 0) error("invalid input argument");

  int perm = 0;
  num_t det = 1;
  for (int i=0, j=0; i < n; i++, j+=n+1)
    det *= a[j], perm += ipiv[i] != i+1;
  mad_free_tmp(a);
  *r = perm & 1 ? -det : det;
  return info;
}

int
mad_cmat_det (const cnum_t x[], cnum_t *r, ssz_t n)
{
  CHKX;
  const int nn=n;
  int info=0, ipiv[n];
  mad_alloc_tmp(cnum_t, a, n*n);
  mad_cvec_copy(x, a, n*n);
  zgetrf_(&nn, &nn, a, &nn, ipiv, &info);

  if (info < 0) error("invalid input argument");

  int perm = 0;
  cnum_t det = 1;
  for (int i=0, j=0; i < n; i++, j+=n+1)
    det *= a[j], perm += ipiv[i] != i+1;
  mad_free_tmp(a);
  *r = perm & 1 ? -det : det;
  return info;
}

// -- inverse -----------------------------------------------------------------o

int
mad_mat_invn (const num_t y[], num_t x, num_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(num_t, u, n*n);
  mad_mat_ident(u, n, n, n);
  int rank = mad_mat_div(u, y, r, n, m, n, rcond);
  mad_free_tmp(u);
  if (x != 1.0) mad_vec_muln(r, x, r, m*n);
  return rank;
}

int // without complex-by-value version
mad_mat_invc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t m, ssz_t n, num_t rcond)
{ CNUM(x); return mad_mat_invc(y, x, r, m, n, rcond); }

int
mad_mat_invc (const num_t y[], cnum_t x, cnum_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(num_t, u, n*n);
  mad_mat_ident(u, n, n, n);
  mad_alloc_tmp(num_t, t, m*n);
  int rank = mad_mat_div(u, y, t, n, m, n, rcond);
  mad_free_tmp(u);
  if (x != 1.0) mad_vec_mulc(t, x, r, m*n);
  mad_free_tmp(t);
  return rank;
}

int
mad_cmat_invn (const cnum_t y[], num_t x, cnum_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(cnum_t, u, n*n);
  mad_cmat_ident(u, n, n, n);
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
  mad_free_tmp(u);
  if (x != 1.0) mad_cvec_muln(r, x, r, m*n);
  return rank;
}

int
mad_cmat_invc (const cnum_t y[], cnum_t x, cnum_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(cnum_t, u, n*n);
  mad_cmat_ident(u, n, n, n);
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
  mad_free_tmp(u);
  if (x != 1.0) mad_cvec_mulc(r, x, r, m*n);
  return rank;
}

int
mad_cmat_invc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t m, ssz_t n, num_t rcond)
{ CNUM(x); return mad_cmat_invc(y, x, r, m, n, rcond); }

// -- divide ------------------------------------------------------------------o

// note:
// X/Y => X * Y^-1 => [m x p] * [p x n] => X:[m x p], Y:[n x p]
// X/Y => X * Y^-1 => (Y'^-1 * X')' => A=Y' and B=X'
// Solving A*X=B => X = A^-1 B = (B'/A')' (col-major!)
//    with Y':[p x n] = A:[M=p x N=n],
//    and  X':[p x m] = B:[M=p x NRHS=m], ipiv:[N]

int
mad_mat_div (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(num_t, a, n*p);
  mad_vec_copy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_vec_copy(x, r, m*p);
    dgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  num_t sz;
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(num_t, rr, ldb*nm);
  mad_mat_copy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  dgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, &info); // compute
  mad_mat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("invalid input argument");
  if (info > 0) error("unexpect lapack error");

  return rank;
}

int
mad_mat_divm (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cnum_t, a, n*p);
  mad_cvec_copy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_vec_copyv(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cnum_t, rr, ldb*nm);
  mad_mat_copym(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("invalid input argument");
  if (info > 0) error("unexpect lapack error");

  return rank;
}

int
mad_cmat_div (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cnum_t, a, n*p);
  mad_cvec_copy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_cvec_copy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cnum_t, rr, ldb*nm);
  mad_cmat_copy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("invalid input argument");
  if (info > 0) error("unexpect lapack error");

  return rank;
}

int
mad_cmat_divm (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cnum_t, a, n*p);
  mad_vec_copyv(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_cvec_copy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cnum_t, rr, ldb*nm);
  mad_cmat_copy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("invalid input argument");
  if (info > 0) error("unexpect lapack error");

  return rank;
}

// -- SVD ---------------------------------------------------------------------o

// SVD decomposition A = U * S * V.t()
// A:[m x n], U:[m x m], S:[min(m,n)], V:[n x n]

int
mad_mat_svd (const num_t x[], num_t u[], num_t s[], num_t v[], ssz_t m, ssz_t n)
{
  assert( x && u && s && v );
  int info=0;
  const int nm=m, nn=n;

  num_t sz;
  int lwork=-1;
  int iwk[8*MIN(m,n)];
  mad_alloc_tmp(num_t, ra, m*n);
  mad_mat_trans(x, ra, m, n);
  dgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn, &sz, &lwork, iwk, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn,  wk, &lwork, iwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra);
  mad_mat_trans(u, u, m, m);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("SVD failed to converged");

  return info;
}

int
mad_cmat_svd (const cnum_t x[], cnum_t u[], num_t s[], cnum_t v[], ssz_t m, ssz_t n)
{
  assert( x && u && s && v );
  int info=0;
  const int nm=m, nn=n;

  cnum_t sz;
  int lwork=-1;
  int iwk[8*MIN(m,n)];
  ssz_t rwk_sz = MIN(m,n) * MAX(5*MIN(m,n)+7, 2*MAX(m,n)+2*MIN(m,n)+1);
  mad_alloc_tmp(num_t, rwk, rwk_sz);
  mad_alloc_tmp(cnum_t, ra, m*n);
  mad_cmat_trans(x, ra, m, n);
  zgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn, &sz, &lwork, rwk, iwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn,  wk, &lwork, rwk, iwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra); mad_free_tmp(rwk);
  mad_cmat_trans(u, u, m, m);
  mad_cvec_conj (v, v, n*n);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("SVD failed to converged");

  return info;
}

// -- EIGEN -------------------------------------------------------------------o

// Eigen values and vectors
// A:[n x n], U:[m x m], S:[min(m,n)], V:[n x n]

int
mad_mat_eigen (const num_t x[], cnum_t w[], num_t vl[], num_t vr[], ssz_t n)
{
  assert( x && w && vl && vr );
  int info=0;
  const int nn=n;

  num_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, wr, n);
  mad_alloc_tmp(num_t, wi, n);
  mad_alloc_tmp(num_t, ra, n*n);
  mad_mat_trans(x, ra, n, n);
  dgeev_("V", "V", &nn, ra, &nn, wr, wi, vl, &nn, vr, &nn, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgeev_("V", "V", &nn, ra, &nn, wr, wi, vl, &nn, vr, &nn,  wk, &lwork, &info); // compute
  mad_vec_cvec(wr, wi, w, n);
  mad_free_tmp(wk); mad_free_tmp(ra);
  mad_free_tmp(wi); mad_free_tmp(wr);
  mad_mat_trans(vl, vl, n, n);
  mad_mat_trans(vr, vr, n, n);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("eigen failed to compute all eigenvalues");

  return info;
}

int
mad_cmat_eigen (const cnum_t x[], cnum_t w[], cnum_t vl[], cnum_t vr[], ssz_t n)
{
  assert( x && w && vl && vr );
  int info=0;
  const int nn=n;

  cnum_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, rwk, 2*n);
  mad_alloc_tmp(cnum_t, ra, n*n);
  mad_cmat_trans(x, ra, n, n);
  zgeev_("V", "V", &nn, ra, &nn, w, vl, &nn, vr, &nn, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgeev_("V", "V", &nn, ra, &nn, w, vl, &nn, vr, &nn,  wk, &lwork, rwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra); mad_free_tmp(rwk);
  mad_cmat_trans(vl, vl, n, n);
  mad_cmat_trans(vr, vr, n, n);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("eigen failed to compute all eigenvalues");

  return info;
}
