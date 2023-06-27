/*
 o-----------------------------------------------------------------------------o
 |
 | Matrix module implementation
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

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"

// --- helpers for debug ------------------------------------------------------o

#if 0
#include <stdio.h>

static void __attribute__((unused))
mprint(str_t name, const num_t a[], ssz_t m, ssz_t n)
{
  printf("%s[%dx%d]=\n", name, m, n);
  FOR(i,m) {
  FOR(j,n) printf("% -10.5f ", a[i*n+j]);
           printf("\n");
  }
}

static void __attribute__((unused))
iprint(str_t name, const idx_t a[], ssz_t m, ssz_t n)
{
  printf("%s[%dx%d]=\n", name, m, n);
  FOR(i,m) {
    FOR(j,n) printf("% 3d ", a[i*n+j]);
             printf("\n");
  }
}
#endif

// --- implementation ---------------------------------------------------------o

#define CHKR   assert( r )
#define CHKX   assert( x )
#define CHKXY  assert( x && y )
#define CHKXR  assert( x && r )
#define CHKYR  assert( y && r )
#define CHKXYR assert( x && y && r )
#define CHKXRX assert( x && r && x != r)

#define CPX(re,im) (* (cpx_t*) & (num_t[2]) { re, im })

// --- matrix, cmatrix, imatrix

struct  matrix { ssz_t nr, nc; num_t data[]; };
struct cmatrix { ssz_t nr, nc; cpx_t data[]; };
struct imatrix { ssz_t nr, nc; idx_t data[]; };

// Note: matrix of zero size are forbidden

void mad_mat_reshape (struct matrix *x, ssz_t m, ssz_t n)
{ CHKX; x->nr = MAX(1,m); x->nc = MAX(1,n); }

void mad_cmat_reshape (struct cmatrix *x, ssz_t m, ssz_t n)
{ CHKX; x->nr = MAX(1,m); x->nc = MAX(1,n); }

void mad_imat_reshape (struct imatrix *x, ssz_t m, ssz_t n)
{ CHKX; x->nr = MAX(1,m); x->nc = MAX(1,n); }

// -----

#if 1
// r[m x n] = x[m x p] * y[p x n]
// naive implementation (more efficient on recent superscalar arch!)
#define MMUL() /* mat * mat */ \
  FOR(i,m*n) r[i] = 0; \
  FOR(i,m) FOR(k,p) FOR(j,n) r[i*n+j] += x[i*p+k] * y[k*n+j];
#else
// [m x n] = [m x p] * [p x n]
// portable vectorized general matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define MMUL() { /* mat * mat */ \
  assert(m>0 && n>0 && p>0); \
  if (n & ~7) { \
    FOR(i,m) FOR(j,0,n-7,8) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = \
      r[i*n+j+4] = r[i*n+j+5] = \
      r[i*n+j+6] = r[i*n+j+7] = 0; \
      FOR(k,p) { \
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
  if (n & 4) { \
    idx_t j = n - (n & 7); \
    FOR(i,m) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = 0; \
      FOR(k,p) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
        r[i*n+j+2] += x[i*p+k] * y[k*n+j+2]; \
        r[i*n+j+3] += x[i*p+k] * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    idx_t j = n - (n & 3); \
    FOR(i,m) { \
      r[i*n+j] = r[i*n+j+1] = 0; \
      FOR(k,p) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    idx_t j = n - 1; \
    FOR(i,m) { \
      r[i*n+j] = 0; \
      FOR(k,p) r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
  } \
}
#endif

// n=1: [m x 1] = [m x p] * [p x 1]
#define MULV() /* mat * vec */ \
  FOR(i,m) r[i] = 0; \
  FOR(i,m) FOR(k,p) r[i] += x[i*p+k] * y[k];

// m=1: [1 x n] = [1 x p] * [p x n]
#define VMUL() /* vec * mat */ \
  FOR(j,n) r[j] = 0; \
  FOR(k,p) FOR(j,n) r[j] += y[k*n+j] * x[k];

// m=1, n=1: [1 x 1] = [1 x p] * [p x 1]
#define IMUL() /* vec * vec */ \
  r[0] = 0; \
  FOR(k,p) r[0] += x[k] * y[k];

// [m x n] = [m x p] * [p x n]
#define MUL() \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: MMUL(); break; \
    case 1: MULV(); break; \
    case 2: VMUL(); break; \
    case 3: IMUL(); break; \
  }

// -----

#if 1
// [m x n] = [p x m]' * [p x n]
// naive implementation (more efficient on recent superscalar arch!)
#define TMMUL(C) /* mat' * mat */ \
  FOR(i,m*n) r[i] = 0; \
  FOR(i,m) FOR(k,p) FOR(j,n) r[i*n+j] += C(x[k*m+i]) * y[k*n+j];
#else
// [m x n] = [p x m]' * [p x n]
// portable vectorized general transpose matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define TMMUL(C) { /* mat' * mat */ \
  assert(m>0 && n>0 && p>0); \
  if (n & ~7) { \
    FOR(i,m) FOR(j,0,n-7,8) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = \
      r[i*n+j+4] = r[i*n+j+5] = \
      r[i*n+j+6] = r[i*n+j+7] = 0; \
      FOR(k,p) { \
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
  if (n & 4) { \
    idx_t j = n - (n & 7); \
    FOR(i,m) { \
      r[i*n+j  ] = r[i*n+j+1] = \
      r[i*n+j+2] = r[i*n+j+3] = 0; \
      FOR(k,p) { \
        r[i*n+j  ] += C(x[k*m+i]) * y[k*n+j  ]; \
        r[i*n+j+1] += C(x[k*m+i]) * y[k*n+j+1]; \
        r[i*n+j+2] += C(x[k*m+i]) * y[k*n+j+2]; \
        r[i*n+j+3] += C(x[k*m+i]) * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    idx_t j = n - (n & 3); \
    FOR(i,m) { \
      r[i*n+j] = r[i*n+j+1] = 0; \
      FOR(k,p) { \
        r[i*n+j  ] += C(x[k*m+i]) * y[k*n+j  ]; \
        r[i*n+j+1] += C(x[k*m+i]) * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    idx_t j = n - 1; \
    FOR(i,m) { \
      r[i*n+j] = 0; \
      FOR(k,p) r[i*n+j] += C(x[k*m+i]) * y[k*n+j]; \
    } \
  } \
}
#endif

// n=1: [m x 1] = [p x m]' * [p x 1]
#define TMULV(C) /* mat' * vec */ \
  FOR(i,m) r[i] = 0; \
  FOR(k,p) FOR(i,m) r[i] += C(x[k*m+i]) * y[k];

// m=1: [1 x n] = [p x 1]' * [p x n]
#define TVMUL(C) /* vec' * mat */ \
  FOR(j,n) r[j] = 0; \
  FOR(k,p) FOR(j,n) r[j] += C(x[k]) * y[k*n+j];

// m=1, n=1: [1 x 1] = [p x 1]' * [p x 1]
#define TIMUL(C) /* vec' * vec */ \
  r[0]= 0; FOR(k,p) r[0] += C(x[k]) * y[k];

// [m x n] = [p x m]' * [p x n]
#define TMUL(C) \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: TMMUL(C); break; \
    case 1: TMULV(C); break; \
    case 2: TVMUL(C); break; \
    case 3: TIMUL(C); break; \
  }

// -----

#if 1
// [m x n] = [m x p] * [n x p]'
// naive implementation (more efficient on recent superscalar arch!)
#define MMULT(C) /* mat * mat' */ \
  FOR(i,m*n) r[i] = 0; \
  FOR(i,m) FOR(j,n) FOR(k,p) r[i*n+j] += x[i*p+k] * C(y[j*p+k]);
#else
// [m x n] = [m x p] * [n x p]'
// portable vectorized general transpose matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
#define MMULT(C) { /* mat * mat' */ \
  assert(m>0 && n>0 && p>0); \
  FOR(i,m*n) r[i] = 0; \
  if (n & ~7) { \
    FOR(i,m) FOR(j,n) FOR(k,0,p-7,8) { \
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
  if (p & 4) { \
    idx_t k = p - (p & 7); \
    FOR(i,m) FOR(j,n) { \
      r[i*n+j] += x[i*p+k  ] * C(y[j*p+k  ]); \
      r[i*n+j] += x[i*p+k+1] * C(y[j*p+k+1]); \
      r[i*n+j] += x[i*p+k+2] * C(y[j*p+k+2]); \
      r[i*n+j] += x[i*p+k+3] * C(y[j*p+k+3]); \
    } \
  } \
  if (p & 2) { \
    idx_t k = p - (p & 3); \
    FOR(i,m) FOR(j,n) { \
      r[i*n+j] += x[i*p+k  ] * C(y[j*p+k  ]); \
      r[i*n+j] += x[i*p+k+1] * C(y[j*p+k+1]); \
    } \
  } \
  if (p & 1) { \
    idx_t k = p - 1; \
    FOR(i,m) FOR(j,n) r[i*n+j] += x[i*p+k] * C(y[j*p+k]); \
  } \
}
#endif

// n=1: [m x 1] = [m x p] * [1 x p]'
#define MULVT(C) /* mat * vec' */ \
  FOR(i,m) r[i] = 0; \
  FOR(i,m) FOR(k,p) r[i] += x[i*p+k] * C(y[k]);

// m=1: [1 x n] = [1 x p]' * [n x p]'
#define VMULT(C) /* vec * mat' */ \
  FOR(j,n) r[j] = 0; \
  FOR(j,n) FOR(k,p) r[j] += x[k] * C(y[j*p+k]);

// m=1, n=1: [1 x 1] = [1 x p] * [1 x p]'
#define IMULT(C) /* vec * vec' */ \
  r[0]= 0; FOR(k,p) r[0] += x[k] * C(y[k]);

// [m x n] = [m x p] * [n x p]'
#define MULT(C) \
  switch(((m == 1) << 1) & (n == 1)) { \
    case 0: MMULT(C); break; \
    case 1: MULVT(C); break; \
    case 2: VMULT(C); break; \
    case 3: IMULT(C); break; \
  }

// -----

// r[m x n] = diag(x[m x p]) * y[p x n]
// naive implementation (more efficient on recent superscalar arch!)
#define DMUL() /* diag(mat) * mat */ \
  if (p == 1) { \
    FOR(i,m) FOR(j,n) r[i*n+j] = x[i] * y[i*n+j]; \
  } else { \
    FOR(i,m*n) r[i] = 0; \
    FOR(i,MIN(m,p)) FOR(j,n) r[i*n+j] = x[i*p+i] * y[i*n+j]; \
  }

// r[m x n] = x[m x p] * diag(y[p x n])
// naive implementation (more efficient on recent superscalar arch!)
#define MULD() /* mat * diag(mat) */ \
  if (p == 1) { \
    FOR(i,m) FOR(j,n) r[i*n+j] = x[i*n+j] * y[j]; \
  } else { \
    FOR(i,m*n) r[i] = 0; \
    FOR(i,m) FOR(j,MIN(n,p)) r[i*n+j] = x[i*p+j] * y[j*n+j]; \
  }

// -----

// [m x n] transpose
#define TRANS(T,C) \
  if (m == 1 || n == 1) { \
    if (x != r || I != C(I)) \
      FOR(i,m*n) r[i] = C(x[i]); \
  } else if ((const void*)x != (const void*)r) { \
    FOR(i,m) FOR(j,n) r[j*m+i] = C(x[i*n+j]); \
  } else if (m == n) { \
    FOR(i,m) FOR(j,i,n) { \
      T t = C(r[j*m+i]); \
      r[j*m+i] = C(r[i*n+j]); \
      r[i*n+j] = t; \
    } \
  } else { \
    mad_alloc_tmp(T, t, m*n); \
    FOR(i,m) FOR(j,n) t[j*m+i] = C(x[i*n+j]); \
    memcpy(r, t, m*n*sizeof(T)); \
    mad_free_tmp(t); \
  };

// -----

// [m x n] copy [+ op]
#define CPY(OP) FOR(i,m) FOR(j,n) r[i*ldr+j] OP##= x[i*ldx+j]

// [m x n] set [+ op]
#define SET(OP) FOR(i,m) FOR(j,n) r[i*ldr+j] OP##= x

// [m x n] sequence [+ op]
#define SEQ(OP) FOR(i,m) FOR(j,n) r[i*ldr+j] OP##= (i*ldr+j)+x

// [m x n] diagonal [+ op]
#define DIAG(OP) FOR(i, MIN(m,n)) r[i*ldr+i] OP##= x

// --- mat

void mad_mat_eye (num_t r[], num_t v, ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; num_t x = 0; SET(); x = v; DIAG(); }

void mad_mat_copy (const num_t x[], num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXRX; CPY(); }

void mad_mat_copym (const num_t x[], cpx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXR; CPY(); }

void mad_mat_trans (const num_t x[], num_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(num_t,); }

void mad_mat_mul (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MUL(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  MUL();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_mulm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { MUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
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

void mad_mat_tmulm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { TMUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
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

void mad_mat_multm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { MULT(conj); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULT(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_dmul (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { DMUL(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  DMUL();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_dmulm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { DMUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  DMUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_muld (const num_t x[], const num_t y[], num_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MULD(); return; }
  mad_alloc_tmp(num_t, r_, m*n);
  num_t *t = r; r = r_;
  MULD();
  mad_vec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_muldm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (y != r) { MULD(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULD();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_mat_center (num_t x[], ssz_t m, ssz_t n, int d)
{ CHKX; num_t mu;
  switch(d) { // 0=vec, 1=row, 2=col, 3=diag
  case 0:
    mu = 0;
    FOR(k,m*n) mu += x[k];
    mu /= m*n;
    FOR(k,m*n) x[k] -= mu;
    break;
  case 1:
    FOR(i,m) { mu = 0;
      FOR(j,n) mu += x[i*n+j];
      mu /= n;
      FOR(j,n) x[i*n+j] -= mu;
    } break;
  case 2:
    FOR(j,n) { mu = 0;
      FOR(i,m) mu += x[i*n+j];
      mu /= m;
      FOR(i,m) x[i*n+j] -= mu;
    } break;
  case 3:
    mu = 0;
    FOR(k,MIN(m,n)) mu += x[k*n+k];
    mu /= MIN(m,n);
    FOR(k,MIN(m,n)) x[k*n+k] -= mu;
    break;
  default: error("invalid direction");
  }
}

void mad_mat_rev (num_t x[], ssz_t m, ssz_t n, int d)
{ CHKX; num_t t;
  switch(d) { // 0=vec, 1=row, 2=col, 3=diag
  case 0: FOR(i,(m*n)/2)        SWAP(x[i    ], x[         m*n-1-i]    , t); break;
  case 1: FOR(i,m  ) FOR(j,n/2) SWAP(x[i*n+j], x[     i *n +n-1-j]    , t); break;
  case 2: FOR(i,m/2) FOR(j,n  ) SWAP(x[i*n+j], x[(m-1-i)*n +    j]    , t); break;
  case 3: FOR(i,MIN(m,n)/2)     SWAP(x[i*n+i], x[(MIN(m,n)-1-i)*(n+1)], t); break;
  default: error("invalid direction");
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
  if (nroll > 0) FOR(i,0,nm,n) {
    mad_vec_copy(x+i+n-nsz, a      ,   nsz); // end of x to a
    mad_vec_copy(x+i      , x+i+nsz, n-nsz); // shift x right
    mad_vec_copy(a        , x+i    ,   nsz); // a to beginning of x
  } else
  if (nroll < 0) FOR(i,0,nm,n) {
    mad_vec_copy(x+i    , a        ,   nsz); // beginning of x to a
    mad_vec_copy(x+i+nsz, x+i      , n-nsz); // shift x left
    mad_vec_copy(a      , x+i+n-nsz,   nsz); // a to end of x
  }
  mad_free_tmp(a);
}

// -- cmat

void mad_cmat_eye (cpx_t r[], cpx_t v, ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; cpx_t x = 0; SET(); x = v; DIAG(); }

void mad_cmat_eye_r (cpx_t r[], num_t v_re, num_t v_im, ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; mad_cmat_eye(r, CPX(v_re,v_im), m, n, ldr); }

void mad_cmat_copy (const cpx_t x[], cpx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXRX; CPY(); }

void mad_cmat_trans (const cpx_t x[], cpx_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(cpx_t,); }

void mad_cmat_ctrans (const cpx_t x[], cpx_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(cpx_t,conj); }

void mad_cmat_mul (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_mulm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { MUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_tmul (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { TMUL(conj); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  TMUL(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_tmulm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { TMUL(conj); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  TMUL(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_mult (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MULT(conj); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULT(conj);
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_multm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { MULT(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULT();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_dmul (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { DMUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  DMUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_dmulm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { DMUL(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  DMUL();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_muld (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r && y != r) { MULD(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULD();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_muldm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p)
{ CHKXYR;
  if (x != r) { MULD(); return; }
  mad_alloc_tmp(cpx_t, r_, m*n);
  cpx_t *t = r; r = r_;
  MULD();
  mad_cvec_copy(r_, t, m*n);
  mad_free_tmp(r_);
}

void mad_cmat_center (cpx_t x[], ssz_t m, ssz_t n, int d)
{ CHKX; cpx_t mu;
  switch(d) { // 0=vec, 1=row, 2=col, 3=diag
  case 0:
    mu = 0;
    FOR(k,m*n) mu += x[k];
    mu /= m*n;
    FOR(k,m*n) x[k] -= mu;
    break;
  case 1:
    FOR(i,m) { mu = 0;
      FOR(j,n) mu += x[i*n+j];
      mu /= n;
      FOR(j,n) x[i*n+j] -= mu;
    } break;
  case 2:
    FOR(j,n) { mu = 0;
      FOR(i,m) mu += x[i*n+j];
      mu /= m;
      FOR(i,m) x[i*n+j] -= mu;
    } break;
  case 3:
    mu = 0;
    FOR(k,MIN(m,n)) mu += x[k*n+k];
    mu /= MIN(m,n);
    FOR(k,MIN(m,n)) x[k*n+k] -= mu;
    break;
  default: error("invalid direction");
  }
}

void mad_cmat_rev (cpx_t x[], ssz_t m, ssz_t n, int d)
{ CHKX; cpx_t t;
  switch(d) { // 0=vec, 1=row, 2=col, 3=diag
  case 0: FOR(i,(m*n)/2)        SWAP(x[i     ], x[         m*n-1-i], t); break;
  case 1: FOR(i,m  ) FOR(j,n/2) SWAP(x[i*n +j], x[     i *n +n-1-j], t); break;
  case 2: FOR(i,m/2) FOR(j,n  ) SWAP(x[i*n +j], x[(m-1-i)*n +    j], t); break;
  case 3: FOR(i,MIN(m,n)/2)     SWAP(x[i*n +i], x[(m-1-i)*n +    i], t); break;
  default: error("invalid direction");
  }
}

void mad_cmat_roll (cpx_t x[], ssz_t m, ssz_t n, int mroll, int nroll)
{ mad_mat_roll((num_t*)x, m, 2*n, mroll, 2*nroll); }

// --- imat

void mad_imat_eye (idx_t r[], idx_t v, ssz_t m, ssz_t n, ssz_t ldr)
{ CHKR; idx_t x = 0; SET(); x = v; DIAG(); }

void mad_imat_copy (const idx_t x[], idx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXRX; CPY(); }

void mad_imat_copym (const idx_t x[], num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)
{ CHKXR; CPY(); }

void mad_imat_trans (const idx_t x[], idx_t r[], ssz_t m, ssz_t n)
{ CHKXR; TRANS(idx_t,); }

void mad_imat_rev (idx_t x[], ssz_t m, ssz_t n, int d)
{ CHKX; idx_t t;
  switch(d) { // 0=vec, 1=row, 2=col, 3=diag
  case 0: FOR(i,(m*n)/2)        SWAP(x[i     ], x[         m*n-1-i], t); break;
  case 1: FOR(i,m  ) FOR(j,n/2) SWAP(x[i*n +j], x[     i *n +n-1-j], t); break;
  case 2: FOR(i,m/2) FOR(j,n  ) SWAP(x[i*n +j], x[(m-1-i)*n +    j], t); break;
  case 3: FOR(i,MIN(m,n)/2)     SWAP(x[i*n +i], x[(m-1-i)*n +    i], t); break;
  default: error("invalid direction");
  }
}

void
mad_imat_roll (idx_t x[], ssz_t m, ssz_t n, int mroll, int nroll)
{ CHKX; mroll %= m; nroll %= n;
  ssz_t nm = n*m, msz = n*abs(mroll), nsz = abs(nroll);
  ssz_t sz = msz > nsz ? msz : nsz;
  mad_alloc_tmp(idx_t, a, sz);
  if (mroll > 0) { // roll rows
    mad_ivec_copy(x+nm-msz, a    ,    msz); // end of x to a
    mad_ivec_copy(x       , x+msz, nm-msz); // shift x down
    mad_ivec_copy(a       , x    ,    msz); // a to beginning of x
  } else
  if (mroll < 0) {
    mad_ivec_copy(x    , a       ,    msz); // beginning of x to a
    mad_ivec_copy(x+msz, x       , nm-msz); // shift x up
    mad_ivec_copy(a    , x+nm-msz,    msz); // a to end of x
  }
  if (nroll > 0) FOR(i,0,nm,n) {
    mad_ivec_copy(x+i+n-nsz, a      ,   nsz); // end of x to a
    mad_ivec_copy(x+i      , x+i+nsz, n-nsz); // shift x right
    mad_ivec_copy(a        , x+i    ,   nsz); // a to beginning of x
  } else
  if (nroll < 0) FOR(i,0,nm,n) {
    mad_ivec_copy(x+i    , a        ,   nsz); // beginning of x to a
    mad_ivec_copy(x+i+nsz, x+i      , n-nsz); // shift x left
    mad_ivec_copy(a      , x+i+n-nsz,   nsz); // a to end of x
  }
  mad_free_tmp(a);
}

// -- Symplectic matrices -----------------------------------------------------o

// M[2n x 2n] accessed as n blocks of [a b ; c d]

#define a_(x,i,j) x[ i   *n+j  ]
#define b_(x,i,j) x[ i   *n+j+1]
#define c_(x,i,j) x[(i+1)*n+j  ]
#define d_(x,i,j) x[(i+1)*n+j+1]

// -- Symplecticity error, compute M' J M - J ---------------------------------o

num_t mad_mat_symperr (const num_t x[], num_t r_[], ssz_t n, num_t *tol_)
{ CHKX; assert(!(n & 1));
  num_t s=0, s0, s1, s2, s3;
  ssz_t nn = n*n;
  mad_alloc_tmp(num_t, r, nn);
  for (idx_t i = 0; i < n-1; i += 2) {
    // i == j
    s1 = -1, s2 = 1;
    for (idx_t k = 0; k < n-1; k += 2) {
      s1 += a_(x,k,i) * d_(x,k,i) - b_(x,k,i) * c_(x,k,i);
      s2 += b_(x,k,i) * c_(x,k,i) - a_(x,k,i) * d_(x,k,i);
    }
    s += s1*s1 + s2*s2;
    b_(r,i,i) = s1, c_(r,i,i) = s2, a_(r,i,i) = d_(r,i,i) = 0;
    // i < j
    for (idx_t j = i+2; j < n-1; j += 2) {
      s0 = s1 = s2 = s3 = 0;
      for (idx_t k = 0; k < n-1; k += 2) {
        s0 += a_(x,k,i) * c_(x,k,j) - a_(x,k,j) * c_(x,k,i);
        s1 += a_(x,k,i) * d_(x,k,j) - b_(x,k,j) * c_(x,k,i);
        s2 += b_(x,k,i) * c_(x,k,j) - a_(x,k,j) * d_(x,k,i);
        s3 += b_(x,k,i) * d_(x,k,j) - b_(x,k,j) * d_(x,k,i);
      }
      s += 2*(s0*s0 + s1*s1 + s2*s2 + s3*s3);
      a_(r,i,j) =  s0, b_(r,i,j) =  s1, c_(r,i,j) =  s2, d_(r,i,j) =  s3;
      a_(r,j,i) = -s0, b_(r,j,i) = -s2, c_(r,j,i) = -s1, d_(r,j,i) = -s3;
    }
  }
  if (tol_) {
    num_t tol = MAX(0,*tol_) ; *tol_ = 1; // is_symp = true
    for (idx_t i = 0; i < nn; i++) if (fabs(r[i]) > tol) {*tol_ = 0; break;}
  }
  if (r_) mad_vec_copy(r, r_, nn);
  mad_free_tmp(r);
  return sqrt(s);
}

num_t mad_cmat_symperr (const cpx_t x[], cpx_t r_[], ssz_t n, num_t *tol_)
{ CHKX; assert(!(n & 1));
  cpx_t s=0, s0, s1, s2, s3;
  ssz_t nn = n*n;
  mad_alloc_tmp(cpx_t, r, nn);
  for (idx_t i = 0; i < n-1; i += 2) {
    // i == j
    s1 = -1, s2 = 1;
    for (idx_t k = 0; k < n-1; k += 2) {
      s1 += conj(a_(x,k,i)) * d_(x,k,i) - b_(x,k,i) * conj(c_(x,k,i));
      s2 += conj(b_(x,k,i)) * c_(x,k,i) - a_(x,k,i) * conj(d_(x,k,i));
    }
    s += s1*s1 + s2*s2;
    b_(r,i,i) = s1, c_(r,i,i) = s2, a_(r,i,i) = d_(r,i,i) = 0;
    // i < j
    for (idx_t j = i+2; j < n-1; j += 2) {
      s0 = s1 = s2 = s3 = 0;
      for (idx_t k = 0; k < n-1; k += 2) {
        s0 += conj(a_(x,k,i)) * c_(x,k,j) - a_(x,k,j) * conj(c_(x,k,i));
        s1 += conj(a_(x,k,i)) * d_(x,k,j) - b_(x,k,j) * conj(c_(x,k,i));
        s2 += conj(b_(x,k,i)) * c_(x,k,j) - a_(x,k,j) * conj(d_(x,k,i));
        s3 += conj(b_(x,k,i)) * d_(x,k,j) - b_(x,k,j) * conj(d_(x,k,i));
      }
      s += 2*(s0*s0 + s1*s1 + s2*s2 + s3*s3);
      a_(r,i,j) =  s0, b_(r,i,j) =  s1, c_(r,i,j) =  s2, d_(r,i,j) =  s3;
      a_(r,j,i) = -s0, b_(r,j,i) = -s2, c_(r,j,i) = -s1, d_(r,j,i) = -s3;
    }
  }
  if (tol_) {
    num_t tol = MAX(0,*tol_) ; *tol_ = 1; // is_symp = true
    for (idx_t i = 0; i < nn; i++) if (cabs(r[i]) > tol) {*tol_ = 0; break;}
  }
  if (r_) mad_cvec_copy(r, r_, nn);
  mad_free_tmp(r);
  return sqrt(cabs(s));
}

// -- Symplectic conjugate, compute \bar{M} = -J M' J -------------------------o

void mad_mat_sympconj (const num_t x[], num_t r[], ssz_t n)
{ CHKXR; assert(!(n & 1));
  num_t t;
  for (idx_t i = 0; i < n-1; i += 2) {     // 2x2 blocks on diagonal
    t = a_(x,i,i),  a_(r,i,i) =  d_(x,i,i),  d_(r,i,i) = t;
    b_(r,i,i) = -b_(x,i,i),  c_(r,i,i) = -c_(x,i,i);

    for (idx_t j = i+2; j < n-1; j += 2) { // 2x2 blocks off diagonal
      t = a_(x,i,j),  a_(r,i,j) =  d_(x,j,i),  d_(r,j,i) =  t;
      t = b_(x,i,j),  b_(r,i,j) = -b_(x,j,i),  b_(r,j,i) = -t;
      t = c_(x,i,j),  c_(r,i,j) = -c_(x,j,i),  c_(r,j,i) = -t;
      t = d_(x,i,j),  d_(r,i,j) =  a_(x,j,i),  a_(r,j,i) =  t;
    }
  }
}

void mad_cmat_sympconj (const cpx_t x[], cpx_t r[], ssz_t n)
{ CHKXR; assert(!(n & 1));
  cpx_t t;
  for (idx_t i = 0; i < n-1; i += 2) {     // 2x2 blocks on diagonal
    t = a_(x,i,i),  a_(r,i,i) =  conj(d_(x,i,i)),  d_(r,i,i) = conj(t);
    b_(r,i,i) = -conj(b_(x,i,i)),  c_(r,i,i) = -conj(c_(x,i,i));

    for (idx_t j = i+2; j < n-1; j += 2) {   // 2x2 blocks off diagonal
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

LAPACK is the default method for solving linear least squares problems for dense
matrices. When the matrix is square and nonsingular the routines dgesv and zgesv
are used otherwise routines dgelsy and zgelsy are used. Alternate methods dgelsd
and zgelsd based on SVD can be used too. The general least squares problems for
dense matrices uses dgglse and zgglse for the first kind, and dggglm and zggglm
for the second kind.

LAPACK is the default method for computing the entire set of singluar values
and singular vectors. For generalized SVD the routines dgesdd and zgesdd are
used.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. For simple eigenvalues the routines dgeev and zgeev are used. For
generalized eigenvalues the routines dggev and zggev are used.

LAPACK C examples using F77 interface can be found at:
https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/
LAPACK F77 examples with data can be found at:
https://github.com/numericalalgorithmsgroup/LAPACK_Examples/tree/master/examples
*/

// -----
// Decompose A = LU with A[m x n] (generalized)
// -----
void dgetrf_ (const int *m, const int *n, num_t A[], const int *lda,
              int *IPIV, int *info);
void zgetrf_ (const int *m, const int *n, cpx_t A[], const int *lda,
              int *IPIV, int *info);

// -----
// Solve A * X = B with A[n x n], B[n x nrhs] and X[n x nrhs]: search min |} b - Ax ||_2 using LU
// -----
void dgesv_ (const int *n, const int *nrhs, num_t A[], const int *lda,
                                 int *IPIV, num_t B[], const int *ldb, int *info);
void zgesv_ (const int *n, const int *nrhs, cpx_t A[], const int *lda,
                                 int *IPIV, cpx_t B[], const int *ldb, int *info);

// -----
// Solve A * X = B with A[m x n], B[m x nrhs] and X[m x nrhs]: search min || b - Ax ||_2 using QR
// -----
void dgelsy_ (const int *m, const int *n, const int *nrhs,
              num_t A[], const int *lda, num_t B[], const int *ldb,
              int jpvt[], const num_t *rcond, int *rank,
              num_t work[], const int lwork[], int *info);
void zgelsy_ (const int *m, const int *n, const int *nrhs,
              cpx_t A[], const int *lda, cpx_t B[], const int *ldb,
              int jpvt[], const num_t *rcond, int *rank,
              cpx_t work[], const int lwork[], num_t rwork[], int *info);

// -----
// Solve A * X = B with A[m x n], B[m x nrhs] and X[m x nrhs]: search min || b - Ax ||_2 using SVD
// -----
void dgelsd_ (const int *m, const int *n, const int *nrhs,
              num_t A[], const int *lda, num_t B[], const int *ldb,
              num_t S[], const num_t *rcond, int *rank,
              num_t work[], int *lwork, int iwork[], int *info);
void zgelsd_ (const int *m, const int *n, const int *nrhs,
              cpx_t A[], const int *lda, cpx_t B[], const int *ldb,
               num_t S[], const num_t *rcond, int *rank,
              cpx_t work[], int *lwork, num_t rwork[], int iwork[], int *info);

// -----
// LS minimization: min_x || c - A*x ||_2 subject to B*x = d using QR
// -----

void dgglse_ (const int *m, const int *n, const int *p,
              num_t A[], const int *lda, num_t B[], const int *ldb,
              num_t C[], num_t D[], num_t X[],
              num_t work[], int *lwork, int *info);
void zgglse_ (const int *m, const int *n, const int *p,
              cpx_t A[], const int *lda, cpx_t B[], const int *ldb,
              cpx_t C[], cpx_t D[], cpx_t X[],
              cpx_t work[], int *lwork, int *info);

// -----
// LS minimization: min_x || y ||_2 subject to A*x + B*y = d using QR
// -----

void dggglm_ (const int *m, const int *n, const int *p,
              num_t A[], const int *lda, num_t B[], const int *ldb,
              num_t D[], num_t X[], num_t Y[],
              num_t work[], int *lwork, int *info);
void zggglm_ (const int *m, const int *n, const int *p,
              cpx_t A[], const int *lda, cpx_t B[], const int *ldb,
              cpx_t D[], cpx_t X[], cpx_t Y[],
              cpx_t work[], int *lwork, int *info);

// -----
// SVD A[m x n]
// -----
void dgesdd_ (str_t jobz, const int *m, const int *n, num_t A[], const int *lda,
              num_t S[], num_t U[], const int *ldu, num_t VT[], const int *ldvt,
              num_t work[], int *lwork, int iwork[], int *info);
void zgesdd_ (str_t jobz, const int *m, const int *n, cpx_t A[], const int *lda,
              num_t S[], cpx_t U[], const int *ldu, cpx_t VT[], const int *ldvt,
              cpx_t work[], int *lwork, num_t rwork[], int iwork[], int *info);

// -----
// Eigen values/vectors A[n x n]
// -----
void dgeev_ (str_t jobvl, str_t jobvr, const int *n, num_t A[], const int *lda,
             num_t WR[], num_t WI[],
             num_t VL[], const int *ldvl, num_t VR[], const int *ldvr,
             num_t work[], int *lwork, int *info);
void zgeev_ (str_t jobvl, str_t jobvr, const int *n, cpx_t A[], const int *lda,
             cpx_t W[], cpx_t VL[], const int *ldvl, cpx_t VR[], const int *ldvr,
             cpx_t work[], int *lwork, num_t rwork[], int *info);

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

  if (info < 0) error("Det: invalid input argument");

  int perm = 0;
  num_t det = 1;
  for (int i=0, j=0; i < n; i++, j+=n+1)
    det *= a[j], perm += ipiv[i] != i+1;
  mad_free_tmp(a);
  *r = perm & 1 ? -det : det;
  return info;
}

int
mad_cmat_det (const cpx_t x[], cpx_t *r, ssz_t n)
{
  CHKX;
  const int nn=n;
  int info=0, ipiv[n];
  mad_alloc_tmp(cpx_t, a, n*n);
  mad_cvec_copy(x, a, n*n);
  zgetrf_(&nn, &nn, a, &nn, ipiv, &info);

  if (info < 0) error("Det: invalid input argument");

  int perm = 0;
  cpx_t det = 1;
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
  mad_mat_eye(u, 1, n, n, n);
#pragma GCC diagnostic push // remove false-positive
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  int rank = mad_mat_div(u, y, r, n, m, n, rcond);
#pragma GCC diagnostic pop
  mad_free_tmp(u);
  if (x != 1) mad_vec_muln(r, x, r, m*n);
  return rank;
}

int // without complex-by-value version
mad_mat_invc_r (const num_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t m, ssz_t n, num_t rcond)
{ return mad_mat_invc(y, CPX(x_re,x_im), r, m, n, rcond); }

int
mad_mat_invc (const num_t y[], cpx_t x, cpx_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(num_t, t, m*n);
  mad_alloc_tmp(num_t, u, n*n);
  mad_mat_eye(u, 1, n, n, n);
  int rank = mad_mat_div(u, y, t, n, m, n, rcond);
  mad_free_tmp(u);
  if (x != 1) mad_vec_mulc(t, x, r, m*n);
  mad_free_tmp(t);
  return rank;
}

int
mad_cmat_invn (const cpx_t y[], num_t x, cpx_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(cpx_t, u, n*n);
  mad_cmat_eye(u, 1, n, n, n);
#pragma GCC diagnostic push // remove false-positive
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
#pragma GCC diagnostic pop
  mad_free_tmp(u);
  if (x != 1) mad_cvec_muln(r, x, r, m*n);
  return rank;
}

int
mad_cmat_invc (const cpx_t y[], cpx_t x, cpx_t r[], ssz_t m, ssz_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  mad_alloc_tmp(cpx_t, u, n*n);
  mad_cmat_eye(u, 1, n, n, n);
#pragma GCC diagnostic push // remove false-positive
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
#pragma GCC diagnostic pop
  mad_free_tmp(u);
  if (x != 1) mad_cvec_mulc(r, x, r, m*n);
  return rank;
}

int
mad_cmat_invc_r (const cpx_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t m, ssz_t n, num_t rcond)
{ return mad_cmat_invc(y, CPX(x_re,x_im), r, m, n, rcond); }

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
    if (info > 0) warn("Div: singular matrix, no solution found");
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

  if (info < 0) error("Div: invalid input argument");
  if (info > 0) error("Div: unexpected lapack error");

  return rank;
}

int
mad_mat_divm (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cpx_t, a, n*p);
  mad_cvec_copy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_vec_copyv(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
    if (info > 0) warn("Div: singular matrix, no solution found");
  }

  // non-square system or singular square system, use QR or LQ factorization
  cpx_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cpx_t, rr, ldb*nm);
  mad_mat_copym(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("Div: invalid input argument");
  if (info > 0) error("Div: unexpected lapack error");

  return rank;
}

int
mad_cmat_div (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cpx_t, a, n*p);
  mad_cvec_copy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_cvec_copy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
    if (info > 0) warn("Div: singular matrix, no solution found");
  }

  // non-square system or singular square system, use QR or LQ factorization
  cpx_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cpx_t, rr, ldb*nm);
  mad_cmat_copy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("Div: invalid input argument");
  if (info > 0) error("Div: unexpected lapack error");

  return rank;
}

int
mad_cmat_divm (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nm=m, nn=n, np=p;
  mad_alloc_tmp(cpx_t, a, n*p);
  mad_vec_copyv(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) {
    int ipiv[n];
    mad_cvec_copy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return mad_free_tmp(a), n;
    if (info > 0) warn("Div: singular matrix, no solution found");
  }

  // non-square system or singular square system, use QR or LQ factorization
  cpx_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  mad_alloc_tmp(cpx_t, rr, ldb*nm);
  mad_cmat_copy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_copy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  mad_free_tmp(wk); mad_free_tmp(rr); mad_free_tmp(a);

  if (info < 0) error("Div: invalid input argument");
  if (info > 0) error("Div: unexpected lapack error");

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
  const int nm=m, nn=n, mn=MIN(m,n);

  num_t sz;
  int lwork=-1;
  int iwk[8*mn];
  mad_alloc_tmp(num_t, ra, m*n);
  mad_mat_trans(x, ra, m, n);
  dgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn, &sz, &lwork, iwk, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn,  wk, &lwork, iwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra);
  mad_mat_trans(u, u, m, m);

  if (info < 0) error("SVD: invalid input argument");
  if (info > 0) warn ("SVD: failed to converge");

  return info;
}

int
mad_cmat_svd (const cpx_t x[], cpx_t u[], num_t s[], cpx_t v[], ssz_t m, ssz_t n)
{
  assert( x && u && s && v );
  int info=0;
  const int nm=m, nn=n, mn=MIN(m,n);

  cpx_t sz;
  int lwork=-1;
  int iwk[8*mn];
  ssz_t rwk_sz = mn * MAX(5*mn+7, 2*MAX(m,n)+2*mn+1);
  mad_alloc_tmp(num_t, rwk, rwk_sz);
  mad_alloc_tmp(cpx_t, ra, m*n);
  mad_cmat_trans(x, ra, m, n);
  zgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn, &sz, &lwork, rwk, iwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgesdd_("A", &nm, &nn, ra, &nm, s, u, &nm, v, &nn,  wk, &lwork, rwk, iwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra); mad_free_tmp(rwk);
  mad_cmat_trans(u, u, m, m);
  mad_cvec_conj (v, v, n*n);

  if (info < 0) error("SVD: invalid input argument");
  if (info > 0) warn ("SVD: failed to converge");

  return info;
}

// -- LEAST SQUARE Solvers ----------------------------------------------------o

int
mad_mat_solve (const num_t a[], const num_t b[], num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  assert( a && b && x );
  int info=0;
  const int nm=m, nn=n, np=p, mn=MAX(m,n);

  num_t sz;
  int lwork=-1, rank;
  int pvt[nn]; memset(pvt, 0, sizeof pvt);
  mad_alloc_tmp(num_t, ta, m*n);
  mad_alloc_tmp(num_t, tb, mn*p); mad_vec_fill(0, tb+m*p, (mn-m)*p);
  mad_vec_copy (b , tb, m*p);
  mad_mat_trans(tb, tb, mn, p);
  mad_mat_trans(a , ta, m , n);
  dgelsy_(&nm, &nn, &np, ta, &nm, tb, &mn, pvt, &rcond, &rank, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgelsy_(&nm, &nn, &np, ta, &nm, tb, &mn, pvt, &rcond, &rank,  wk, &lwork, &info); // compute
  mad_mat_trans(tb, tb, p, mn);
  mad_vec_copy (tb,  x, n*p);

  mad_free_tmp(wk); mad_free_tmp(ta); mad_free_tmp(tb);

  if (info < 0) error("Solve: invalid input argument");
  if (info > 0) warn ("Solve: unexpected lapack error");

  return rank;
}

int
mad_cmat_solve (const cpx_t a[], const cpx_t b[], cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)
{
  assert( a && b && x );
  int info=0;
  const int nm=m, nn=n, np=p, mn=MAX(m,n);

  cpx_t sz;
  num_t rwk[2*nn];
  int lwork=-1, rank;
  int pvt[nn]; memset(pvt, 0, sizeof pvt);
  mad_alloc_tmp(cpx_t, ta, m*n);
  mad_alloc_tmp(cpx_t, tb, mn*p); mad_cvec_fill(0, tb+m*p, (mn-m)*p);
  mad_cvec_copy (b , tb, m*p);
  mad_cmat_trans(tb, tb, mn, p);
  mad_cmat_trans(a , ta, m , n);
  zgelsy_(&nm, &nn, &np, ta, &nm, tb, &mn, pvt, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgelsy_(&nm, &nn, &np, ta, &nm, tb, &mn, pvt, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_trans(tb, tb, p, mn);
  mad_cvec_copy (tb,  x, n*p);

  mad_free_tmp(wk); mad_free_tmp(ta); mad_free_tmp(tb);

  if (info < 0) error("Solve: invalid input argument");
  if (info > 0) warn ("Solve: unexpected lapack error");

  return rank;
}

int
mad_mat_ssolve (const num_t a[], const num_t b[], num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[])
{
  assert( a && b && x );
  int info=0;
  const int nm=m, nn=n, np=p, mn=MAX(m,n);

  num_t sz;
  int lwork=-1, rank, isz;
  mad_alloc_tmp(num_t, ta, m *n);
  mad_alloc_tmp(num_t, tb, mn*p);
  mad_alloc_tmp(num_t, ts, MIN(m,n));
  mad_vec_copy (b , tb, m*p);
  mad_vec_fill (0 , tb +m*p, (mn-m)*p);
  mad_mat_trans(tb, tb, mn, p);
  mad_mat_trans(a , ta, m , n);
  dgelsd_(&nm, &nn, &np, ta, &nm, tb, &mn, ts, &rcond, &rank, &sz, &lwork, &isz, &info); // query
  mad_alloc_tmp(num_t,  wk, lwork=sz);
  mad_alloc_tmp(int  , iwk, isz);
  dgelsd_(&nm, &nn, &np, ta, &nm, tb, &mn, ts, &rcond, &rank,  wk, &lwork,  iwk, &info); // compute
  mad_mat_trans(tb, tb, p, mn);
  mad_vec_copy (tb,  x, n*p);

  if (s_) mad_vec_copy(ts, s_, MIN(m,n));

  mad_free_tmp(wk); mad_free_tmp(iwk);
  mad_free_tmp(ta); mad_free_tmp(tb); mad_free_tmp(ts);

  if (info < 0) error("SSolve: invalid input argument");
  if (info > 0) warn ("SSolve: failed to converge");

  return rank;
}

int
mad_cmat_ssolve (const cpx_t a[], const cpx_t b[], cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[])
{
  assert( a && b && x );
  int info=0;
  const int nm=m, nn=n, np=p, mn=MAX(m,n);

  num_t rsz;
  cpx_t sz;
  int lwork=-1, rank, isz;
  mad_alloc_tmp(cpx_t, ta, m*n);
  mad_alloc_tmp(cpx_t, tb, mn*p);
  mad_alloc_tmp(num_t, ts, MIN(m,n));
  mad_cvec_copy (b , tb, m*p);
  mad_cvec_fill (0 , tb +m*p, (mn-m)*p);
  mad_cmat_trans(tb, tb, mn, p);
  mad_cmat_trans(a , ta, m , n);
  zgelsd_(&nm, &nn, &np, ta, &nm, tb, &mn, ts, &rcond, &rank, &sz, &lwork, &rsz, &isz, &info); // query
  mad_alloc_tmp(cpx_t,  wk, lwork=creal(sz));
  mad_alloc_tmp( num_t, rwk, (int)rsz);
  mad_alloc_tmp( int  , iwk, isz);
  zgelsd_(&nm, &nn, &np, ta, &nm, tb, &mn, ts, &rcond, &rank,  wk, &lwork,  rwk,  iwk, &info); // compute
  mad_cmat_trans(tb, tb, p, mn);
  mad_cvec_copy (tb,  x, n*p);

  if (s_) mad_vec_copy(ts, s_, MIN(m,n));

  mad_free_tmp(wk); mad_free_tmp(rwk); mad_free_tmp(iwk);
  mad_free_tmp(ta); mad_free_tmp(tb);  mad_free_tmp(ts);

  if (info < 0) error("SSolve: invalid input argument");
  if (info > 0) warn ("SSolve: failed to converge");

  return rank;
}

// -- Generalized LS Solvers --------------------------------------------------o

int
mad_mat_gsolve (const num_t a[], const num_t b[], const num_t c[], const num_t d[],
                num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_)
{
  assert( a && b && x );
  ensure( 0 <= p && p <= n && n <= m+p, "invalid system sizes" );
  int info=0;
  const int nm=m, nn=n, np=p;

  num_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, ta, m*n);
  mad_alloc_tmp(num_t, tb, p*n);
  mad_alloc_tmp(num_t, tc, m);
  mad_alloc_tmp(num_t, td, p);
  mad_mat_trans(a, ta, m, n);
  mad_mat_trans(b, tb, p, n);
  mad_vec_copy (c, tc, m);
  mad_vec_copy (d, td, p);
  dgglse_(&nm, &nn, &np, ta, &nm, tb, &np, tc, td, x, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgglse_(&nm, &nn, &np, ta, &nm, tb, &np, tc, td, x,  wk, &lwork, &info); // compute

  if (nrm_) *nrm_ = mad_vec_norm(tc+(n-p), m-(n-p)); // residues

  mad_free_tmp(wk);
  mad_free_tmp(ta); mad_free_tmp(tb); mad_free_tmp(tc); mad_free_tmp(td);

  if (info < 0) error("GSolve: invalid input argument");
  if (info > 0) warn ("GSolve: [B A] is singular, no solution found");

  return info;
}

int
mad_cmat_gsolve (const cpx_t a[], const cpx_t b[], const cpx_t c[], const cpx_t d[],
                 cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_)
{
  assert( a && b && x );
  ensure( 0 <= p && p <= n && n <= m+p, "invalid system sizes" );
  int info=0;
  const int nm=m, nn=n, np=p;

  cpx_t sz;
  int lwork=-1;
  mad_alloc_tmp(cpx_t, ta, m*n);
  mad_alloc_tmp(cpx_t, tb, p*n);
  mad_alloc_tmp(cpx_t, tc, m);
  mad_alloc_tmp(cpx_t, td, p);
  mad_cmat_trans(a, ta, m, n);
  mad_cmat_trans(b, tb, p, n);
  mad_cvec_copy (c, tc, m);
  mad_cvec_copy (d, td, p);
  zgglse_(&nm, &nn, &np, ta, &nm, tb, &np, tc, td, x, &sz, &lwork, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=sz);
  zgglse_(&nm, &nn, &np, ta, &nm, tb, &np, tc, td, x,  wk, &lwork, &info); // compute

  if (nrm_) *nrm_ = mad_cvec_norm(tc+(n-p), m-(n-p)); // residues

  mad_free_tmp(wk);
  mad_free_tmp(ta); mad_free_tmp(tb); mad_free_tmp(tc); mad_free_tmp(td);

  if (info < 0) error("GSolve: invalid input argument");
  if (info > 0) warn ("GSolve: [B A] is singular, no solution found");

  return info;
}

int
mad_mat_gmsolve (const num_t a[], const num_t b[], const num_t d[],
                 num_t x[], num_t y[], ssz_t m, ssz_t n, ssz_t p)
{
  assert( a && b && x );
  ensure( 0 <= p && n <= m && m <= n+p, "invalid system sizes" );
  int info=0;
  const int nm=m, nn=n, np=p;

  num_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, ta, m*n);
  mad_alloc_tmp(num_t, tb, m*p);
  mad_alloc_tmp(num_t, td, m);
  mad_mat_trans(a, ta, m, n);
  mad_mat_trans(b, tb, m, p);
  mad_vec_copy (d, td, m);
  dggglm_(&nm, &nn, &np, ta, &nm, tb, &nm, td, x, y, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dggglm_(&nm, &nn, &np, ta, &nm, tb, &nm, td, x, y,  wk, &lwork, &info); // compute

  mad_free_tmp(wk);
  mad_free_tmp(ta); mad_free_tmp(tb); mad_free_tmp(td);

  if (info < 0) error("GMSolve: invalid input argument");
  if (info > 0) warn ("GMSolve: [A B] is singular, no solution found");

  return info;
}

int
mad_cmat_gmsolve (const cpx_t a[], const cpx_t b[], const cpx_t d[],
                  cpx_t x[], cpx_t y[], ssz_t m, ssz_t n, ssz_t p)
{
  assert( a && b && x );
  ensure( 0 <= p && n <= m && m <= n+p, "invalid system sizes" );
  int info=0;
  const int nm=m, nn=n, np=p;

  cpx_t sz;
  int lwork=-1;
  mad_alloc_tmp(cpx_t, ta, m*n);
  mad_alloc_tmp(cpx_t, tb, m*p);
  mad_alloc_tmp(cpx_t, td, m);
  mad_cmat_trans(a, ta, m, n);
  mad_cmat_trans(b, tb, m, p);
  mad_cvec_copy (d, td, m);
  zggglm_(&nm, &nn, &np, ta, &nm, tb, &nm, td, x, y, &sz, &lwork, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=sz);
  zggglm_(&nm, &nn, &np, ta, &nm, tb, &nm, td, x, y,  wk, &lwork, &info); // compute

  mad_free_tmp(wk);
  mad_free_tmp(ta); mad_free_tmp(tb); mad_free_tmp(td);

  if (info < 0) error("GMSolve: invalid input argument");
  if (info > 0) warn ("GMSolve: [A B] is singular, no solution found");

  return info;
}

// -- EIGEN -------------------------------------------------------------------o

// Eigen values and vectors
// A:[n x n], U:[m x m], S:[min(m,n)], V:[n x n]

int
mad_mat_eigen (const num_t x[], cpx_t w[], num_t vl_[], num_t vr_[], ssz_t n)
{
  assert( x && w );
  int info=0;
  const int nn=n;
  const str_t vls = vl_ ? "V" : "N";
  const str_t vrs = vr_ ? "V" : "N";

  num_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, wr, n);
  mad_alloc_tmp(num_t, wi, n);
  mad_alloc_tmp(num_t, ra, n*n);
  mad_mat_trans(x, ra, n, n);
  dgeev_(vls, vrs, &nn, ra, &nn, wr, wi, vl_, &nn, vr_, &nn, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dgeev_(vls, vrs, &nn, ra, &nn, wr, wi, vl_, &nn, vr_, &nn,  wk, &lwork, &info); // compute
  mad_vec_cplx(wr, wi, w, n);
  mad_free_tmp(wk); mad_free_tmp(ra);
  mad_free_tmp(wi); mad_free_tmp(wr);
//if (vl_) mad_mat_trans(vl_, vl_, n, n);
  if (vr_) mad_mat_trans(vr_, vr_, n, n);

  if (info < 0) error("Eigen: invalid input argument");
  if (info > 0) warn ("Eigen: failed to compute all eigenvalues");

  return info;
}

int
mad_cmat_eigen (const cpx_t x[], cpx_t w[], cpx_t vl_[], cpx_t vr_[], ssz_t n)
{
  assert( x && w );
  int info=0;
  const int nn=n;
  const str_t vls = vl_ ? "V" : "N";
  const str_t vrs = vr_ ? "V" : "N";

  cpx_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, rwk, 2*n);
  mad_alloc_tmp(cpx_t, ra, n*n);
  mad_cmat_trans(x, ra, n, n);
  zgeev_(vls, vrs, &nn, ra, &nn, w, vl_, &nn, vr_, &nn, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cpx_t, wk, lwork=creal(sz));
  zgeev_(vls, vrs, &nn, ra, &nn, w, vl_, &nn, vr_, &nn,  wk, &lwork, rwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra); mad_free_tmp(rwk);
//if (vl_) mad_cmat_trans(vl_, vl_, n, n);
  if (vr_) mad_cmat_trans(vr_, vr_, n, n);

  if (info < 0) error("Eigen: invalid input argument");
  if (info > 0) warn ("Eigen: failed to compute all eigenvalues");

  return info;
}

// -- GEOMETRY ----------------------------------------------------------------o

// -- Helpers -----------------------------------------------------------------o

#define NN            (N*N)
#define X(i,j)        x[(i-1)*N+(j-1)]
#define VCPY(src,dst) for(idx_t i=0; i<N ; dst[i]=src[i], ++i)
#define MCPY(src,dst) for(idx_t i=0; i<NN; dst[i]=src[i], ++i)

// -- 2D geometry -------------------------------------------------------------o

#define N 2

// 2D rotation

void mad_mat_rot (num_t x[NN], num_t a) // R
{
  CHKX;
  num_t ca = cos(a), sa = sin(a);
  num_t r[NN] = {ca,-sa,
                 sa, ca};
  MCPY(r,x);
}

#undef N

// -- 3D geometry -------------------------------------------------------------o

#define N 3

// 3D rotations (one axis)

void mad_mat_rotx (num_t x[NN], num_t ax) // Rx
{
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t r[NN] = {1,  0,  0,
                 0, cx,-sx,
                 0, sx, cx};
  MCPY(r,x);
}

void mad_mat_roty (num_t x[NN], num_t ay) // Ry
{
  CHKX;
  num_t cy = cos(ay), sy = sin(ay);
  num_t r[NN] = { cy, 0, sy,
                   0, 1,  0,
                 -sy, 0, cy};
  MCPY(r,x);
}

void mad_mat_rotz (num_t x[NN], num_t az) // Rz
{
  CHKX;
  num_t cz = cos(az), sz = sin(az);
  num_t r[NN] = {cz,-sz, 0,
                 sz, cz, 0,
                  0,  0, 1};
  MCPY(r,x);
}

// 3D rotations (two axis)

void mad_mat_rotxy (num_t x[NN], num_t ax, num_t ay, log_t inv) // Ry.Rx
{
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t cy = cos(ay), sy = sin(ay);

  if (!inv) {  // normal
    num_t r[NN] = { cy, sx*sy, cx*sy,
                     0,    cx,   -sx,
                   -sy, sx*cy, cx*cy};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {   cy,   0,   -sy,
                   sx*sy,  cx, sx*cy,
                   cx*sy, -sx, cx*cy};
    MCPY(r,x);
  }
}

void mad_mat_rotxz (num_t x[NN], num_t ax, num_t az, log_t inv) // Rz.Rx
{
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t cz = cos(az), sz = sin(az);

  if (!inv) {  // normal
    num_t r[NN] = {cz,-cx*sz, sx*sz,
                   sz, cx*cz,-sx*cz,
                    0,    sx,    cx};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {    cz,    sz,  0,
                   -cx*sz, cx*cz, sx,
                    sx*sz,-sx*cz, cx};
    MCPY(r,x);
  }
}

void mad_mat_rotyz (num_t x[NN], num_t ay, num_t az, log_t inv) // Rz.Ry
{
  CHKX;
  num_t cy = cos(ay), sy = sin(ay);
  num_t cz = cos(az), sz = sin(az);

  if (!inv) {  // normal
    num_t r[NN] = {cy*cz,-sz, sy*cz,
                   cy*sz, cz, sy*sz,
                     -sy,  0,    cy};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {cy*cz,cy*sz, -sy,
                     -sz,   cz,   0,
                   sy*cz,sy*sz,  cy};
    MCPY(r,x);
  }
}

// 3D rotations (three axis)

void mad_mat_rotxyz (num_t x[NN], num_t ax, num_t ay, num_t az, log_t inv)
{ // Rz.Ry.Rx
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t cy = cos(ay), sy = sin(ay);
  num_t cz = cos(az), sz = sin(az);

  if (!inv) {  // normal
    num_t r[NN] = {cy*cz, cz*sx*sy - cx*sz, cx*cz*sy + sx*sz,
                   cy*sz, sx*sy*sz + cx*cz, cx*sy*sz - cz*sx,
                     -sy,            cy*sx,            cx*cy};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {           cy*cz,            cy*sz,   -sy,
                   cz*sx*sy - cx*sz, sx*sy*sz + cx*cz, cy*sx,
                   cx*cz*sy + sx*sz, cx*sy*sz - cz*sx, cx*cy};
    MCPY(r,x);
  }
}

void mad_mat_rotxzy (num_t x[NN], num_t ax, num_t ay, num_t az, log_t inv)
{ // Ry.Rz.Rx
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t cy = cos(ay), sy = sin(ay);
  num_t cz = cos(az), sz = sin(az);

  if (!inv) {  // normal
    num_t r[NN] = { cy*cz, sx*sy - cx*cy*sz,  cx*sy + cy*sx*sz,
                       sz, cx*cz           , -cz*sx           ,
                   -cz*sy, cy*sx + cx*sy*sz,  cx*cy - sx*sy*sz};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {cy*cz           ,     sz, -cz*sy           ,
                   sx*sy - cx*cy*sz,  cx*cz,  cy*sx + cx*sy*sz,
                   cx*sy + cy*sx*sz, -cz*sx,  cx*cy - sx*sy*sz};
    MCPY(r,x);
  }
}

void mad_mat_rotyxz (num_t x[NN], num_t ax, num_t ay, num_t az, log_t inv)
{ // Rz.Rx.Ry
  CHKX;
  num_t cx = cos(ax), sx = sin(ax);
  num_t cy = cos(ay), sy = sin(ay);
  num_t cz = cos(az), sz = sin(az);

  if (!inv) {  // normal
    num_t r[NN] = { cy*cz - sx*sy*sz, -cx*sz, cz*sy + cy*sx*sz,
                    cy*sz + cz*sx*sy,  cx*cz, sy*sz - cy*cz*sx,
                   -cx*sy           ,     sx, cx*cy           };
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = { cy*cz - sx*sy*sz, cy*sz + cz*sx*sy, -cx*sy,
                   -cx*sz           , cx*cz           ,     sx,
                    cz*sy + cy*sx*sz, sy*sz - cy*cz*sx,  cx*cy};
    MCPY(r,x);
  }
}

// 3D angles from rotations

void mad_mat_torotxyz (const num_t x[NN], num_t r[N], log_t inv)
{ // extract ax, ay, az from rotxyz
  CHKXR;

  num_t x11 = X(1,1), x33 = X(3,3), x21, x31, x32;

  if (!inv) x21 = X(2,1), x31 = X(3,1), x32 = X(3,2);
  else      x21 = X(1,2), x31 = X(1,3), x32 = X(2,3);

  r[0] = atan2( x32, x33 );                     // ax
  r[1] = atan2(-x31, sqrt(x32*x32 + x33*x33) ); // ay
  r[2] = atan2( x21, x11 );                     // az
}

void mad_mat_torotxzy (const num_t x[NN], num_t r[N], log_t inv)
{ // extract ax, ay, az from rotxzy
  CHKXR;
  num_t x11 = X(1,1), x22 = X(2,2), x21, x23, x31;

  if (!inv) x21 = X(2,1), x23 = X(2,3), x31 = X(3,1);
  else      x21 = X(1,2), x23 = X(3,2), x31 = X(1,3);

  r[0] = atan2(-x23, x22 );                     // ax
  r[1] = atan2(-x31, x11 );                     // ay
  r[2] = atan2( x21, sqrt(x22*x22 + x23*x23) ); // az
}

void mad_mat_torotyxz (const num_t x[NN], num_t r[N], log_t inv)
{ // extract ax, ay, az from rotyxz
  CHKXR;
  num_t x22 = X(2,2), x33 = X(3,3), x12, x31, x32;

  if (!inv) x12 = X(1,2), x31 = X(3,1), x32 = X(3,2);
  else      x12 = X(2,1), x31 = X(1,3), x32 = X(2,3);

  r[0] = atan2( x32, sqrt(x12*x12 + x22*x22) ); // ax
  r[1] = atan2(-x31, x33 );                     // ay
  r[2] = atan2(-x12, x22 );                     // az
}

// 3D vector rotation

void mad_mat_rotv (num_t x[NN], const num_t v[N], num_t a, log_t inv)
{
  assert(x && v);

  num_t vx = v[0], vy = v[1], vz = v[2];
  num_t n = vx*vx + vy*vy + vz*vz;

  if (n == 0) {
    mad_mat_eye(x, 1, N, N, N);
    return;
  }

  if (n != 1) {
    n = 1/sqrt(n);
    vx *= n, vy *= n, vz *= n;
  }

  num_t xx = vx*vx,  yy = vy*vy,  zz = vz*vz;
  num_t xy = vx*vy,  xz = vx*vz,  yz = vy*vz;
  num_t ca = cos(a), sa = sin(a), C  = 1-ca;

  if (!inv) {  // normal
    num_t r[NN] = {xx*C +    ca, xy*C - vz*sa, xz*C + vy*sa,
                   xy*C + vz*sa, yy*C +    ca, yz*C - vx*sa,
                   xz*C - vy*sa, yz*C + vx*sa, zz*C +    ca};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {xx*C +    ca, xy*C + vz*sa, xz*C - vy*sa,
                   xy*C - vz*sa, yy*C +    ca, yz*C + vx*sa,
                   xz*C + vy*sa, yz*C - vx*sa, zz*C +    ca};
    MCPY(r,x);
  }
}

num_t mad_mat_torotv (const num_t x[NN], num_t v_[N], log_t inv)
{
  CHKX;
  num_t vx, vy, vz;

  if (!inv) {
    vx = X(3,2) - X(2,3);
    vy = X(1,3) - X(3,1);
    vz = X(2,1) - X(1,2);
  } else {
    vx = X(2,3) - X(3,2);
    vy = X(3,1) - X(1,3);
    vz = X(1,2) - X(2,1);
  }

  num_t n = sqrt(vx*vx + vy*vy + vz*vz);
  num_t t = X(1,1) + X(2,2) + X(3,3);
  num_t a = atan2(n, t-1);

  if (v_) {
    n = n != 0 ? 1/n : 0;
    v_[0] = n*vx, v_[1] = n*vy, v_[2] = n*vz;
  }
  return a;
}

// Quaternion

void mad_mat_rotq (num_t x[NN], const num_t q[4], log_t inv)
{
  assert(x && q);

  num_t qw = q[0], qx = q[1], qy = q[2], qz = q[3];
  num_t n = qw*qw + qx*qx + qy*qy + qz*qz;
  num_t s = n != 0 ? 2/n : 0;
  num_t wx = s*qw*qx, wy = s*qw*qy, wz = s*qw*qz;
  num_t xx = s*qx*qx, xy = s*qx*qy, xz = s*qx*qz;
  num_t yy = s*qy*qy, yz = s*qy*qz, zz = s*qz*qz;

  if (!inv) {  // normal
    num_t r[NN] = {1-(yy+zz),    xy-wz ,    xz+wy,
                      xy+wz , 1-(xx+zz),    yz-wx,
                      xz-wy ,    yz+wx , 1-(xx+yy)};
    MCPY(r,x);
  } else {     // transposed
    num_t r[NN] = {1-(yy+zz),    xy+wz ,    xz-wy,
                      xy-wz , 1-(xx+zz),    yz+wx,
                      xz+wy ,    yz-wx , 1-(xx+yy)};
    MCPY(r,x);
  }
}

void mad_mat_torotq (const num_t x[NN], num_t q[4], log_t inv)
{
  CHKX;
  num_t xx = X(1,1), yy = X(2,2), zz = X(3,3);
  num_t tt = xx+yy+zz, rr, ss;

  // stable trace
  if (tt > -0.99999) {
    rr = sqrt(1+tt), ss = 0.5/rr;
    q[0] = 0.5*rr;
    if (!inv) {  // normal
      q[1] = (X(3,2) - X(2,3)) * ss;
      q[2] = (X(1,3) - X(3,1)) * ss;
      q[3] = (X(2,1) - X(1,2)) * ss;
    } else {     // transposed
      q[1] = (X(2,3) - X(3,2)) * ss;
      q[2] = (X(3,1) - X(1,3)) * ss;
      q[3] = (X(1,2) - X(2,1)) * ss;
    }
    return;
  }

  // look for more stable trace
  num_t m = MAX(xx, yy, zz);
  if (!inv) {  // normal
    if (m == xx) {
      rr = sqrt(1+xx-yy-zz), ss = 0.5/rr;
      q[1] = 0.5*rr;
      q[0] = (X(3,2) - X(2,3)) * ss;
      q[2] = (X(1,3) + X(3,1)) * ss;
      q[3] = (X(2,1) + X(1,2)) * ss;
    } else if (m == yy) {
      rr = sqrt(1+yy-xx-zz), ss = 0.5/rr;
      q[2] = 0.5*rr;
      q[0] = (X(3,2) - X(2,3)) * ss;
      q[1] = (X(1,3) - X(3,1)) * ss;
      q[3] = (X(2,1) + X(1,2)) * ss;
    } else {
      rr = sqrt(1+zz-xx-yy), ss = 0.5/rr;
      q[3] = 0.5*rr;
      q[0] = (X(3,2) - X(2,3)) * ss;
      q[1] = (X(1,3) - X(3,1)) * ss;
      q[2] = (X(2,1) - X(1,2)) * ss;
    }
  } else {     // transposed
    if (m == xx) {
      rr = sqrt(1+xx-yy-zz), ss = 0.5/rr;
      q[1] = 0.5*rr;
      q[0] = (X(2,3) - X(3,2)) * ss;
      q[2] = (X(3,1) + X(1,3)) * ss;
      q[3] = (X(1,2) + X(2,1)) * ss;
    } else if (m == yy) {
      rr = sqrt(1+yy-xx-zz), ss = 0.5/rr;
      q[2] = 0.5*rr;
      q[0] = (X(2,3) - X(3,2)) * ss;
      q[1] = (X(3,1) - X(1,3)) * ss;
      q[3] = (X(1,2) + X(2,1)) * ss;
    } else {
      rr = sqrt(1+zz-xx-yy), ss = 0.5/rr;
      q[3] = 0.5*rr;
      q[0] = (X(2,3) - X(3,2)) * ss;
      q[1] = (X(3,1) - X(1,3)) * ss;
      q[2] = (X(1,2) - X(2,1)) * ss;
    }
  }
}

#undef N
#undef X

// -- Orbit Correction --------------------------------------------------------o

int mad_use_madx_micado = 0;
int mad_use_madx_svdcnd = 0;

extern void // see madx_micado.f90
micit_(num_t cin[], num_t res[],
       int nx[], num_t *rms, int *im, int *ic, int *iter,
       /* working buffers */
       int ny[], num_t ax[], num_t cinx[], num_t xinx[], num_t resx[],
       num_t rho[], num_t ptop[], num_t rmss[], num_t xrms[], num_t xptp[],
       num_t xiter[], int *ifail);

extern void
svddec_(num_t svdmat[], num_t umat[], num_t vmat[],
        /* working buffers */
        num_t ws[], num_t wvec[], int sortw[], num_t *sngcut, num_t *sngval,
        /* sizes and output */
        int *im, int *ic, int *iflag, int sing[]);

static void
vec_sort (num_t v[], idx_t c[], ssz_t n)
{
  num_t tr;
  idx_t ti;

  // Set indexes.
  FOR(i,n) c[i] = i;

  // Sort values by ascending order.
  FOR(i,1,n) RFOR(j,i+1)
    if (v[j-1] > v[j]) {
      SWAP(v[j-1], v[j], tr);
      SWAP(c[j-1], c[j], ti);
    }
}

static ssz_t
ivec_sort (idx_t v[], ssz_t n, log_t rmdup)
{
  idx_t ti;

  // Sort indexes by ascending order.
  FOR(i,1,n) RFOR(j,i+1)
    if (v[j-1] > v[j]) SWAP(v[j-1], v[j], ti);

  // Remove duplicates.
  if (rmdup) {
    idx_t k = 1;
    FOR(i,1,n) if (v[k-1] < v[i]) v[k++] = v[i];
    return k;
  }

  return n;
}

static int // madx legacy code wrapper
madx_svdcnd (const num_t a[], idx_t c[], ssz_t m, ssz_t n, num_t scut, num_t s_[], num_t sval)
{
  /* copy buffers */
  mad_alloc_tmp(num_t, A  , m*n);
  mad_alloc_tmp(num_t, U  , m*n);
  mad_alloc_tmp(num_t, V  , n*n);
  mad_alloc_tmp(num_t, S  , n  );
  /* working buffers */
  mad_alloc_tmp(num_t, W  , n  );
  mad_alloc_tmp(idx_t, srt, n  );
  mad_alloc_tmp(idx_t, sng, 2*n);

  mad_mat_trans(a, A, m, n);

  int im=m, ic=n, nc=0;

  svddec_(A, U, V, W, S, srt, &scut, &sval, &im, &ic, &nc, sng);

  // Backup singular values.
  if (s_) mad_vec_copy(S, s_, MIN(m,n));

  // Backup indexes of columns to remove.
  FOR(i,nc) c[i] = sng[2*i];

  /* copy buffers */
  mad_free_tmp(A);
  mad_free_tmp(U);
  mad_free_tmp(V);
  mad_free_tmp(S);
  /* working buffers */
  mad_free_tmp(W);
  mad_free_tmp(srt);
  mad_free_tmp(sng);

  // Return sorted indexes of columns to remove.
  return ivec_sort(c, nc, true);
}

int // Matrix preconditionning using SVD, return indexes of columns to remove.
mad_mat_svdcnd(const num_t a[], idx_t c[], ssz_t m, ssz_t n,
               ssz_t N, num_t rcond, num_t s_[], num_t tol)
{
  assert(a && c);

  // X-check with MAD-X SVD cond (where N = n-5)
  if (mad_use_madx_svdcnd) {
    return madx_svdcnd(a, c, m, n, rcond, s_, 1/tol);
  }

  ssz_t mn = MIN(m,n);

  mad_alloc_tmp(num_t, U, m*m);
  mad_alloc_tmp(num_t, V, n*n);
  mad_alloc_tmp(num_t, S, mn );

  int info = mad_mat_svd(a, U, S, V, m, n);
  if (info != 0) return -1;

  // Backup singular values.
  if (s_) mad_vec_copy(S, s_, mn);

  // N == 0 means to check for all singular values.
  if (N > mn || N <= 0) N = mn;

  // Tolerance on components similarity in V columns.
  tol = MAX(tol, DBL_EPSILON);

  // Tolerance on keeping singular values.
  rcond = MAX(rcond, DBL_EPSILON);

  // Number of columns to remove.
  idx_t nc = 0;

#define V(i,j) V[(i)*n+(j)]

  // Loop over increasing singular values.
  RFOR(i,mn,mn-N) {
    // Singular value is large, stop checking.
    if (S[i] > rcond*S[0]) break;

    // Loop over rows of V (i.e. columns of V^T)
    FOR(j,n-1) FOR(k,j+1,n) {
      num_t vj = fabs(V(j,i));

      // Proceed only significant component for this singular value.
      if (vj > 1e-4) {
        num_t vk  = fabs(V(k,i));
        num_t rat = fabs(vj-vk)/(vj+vk);

        // Discard column j with similar (or opposite) effect of column k > j.
        if (rat <= tol) {
          c[nc++] = j; // can hold duplicated indexes...
          if (nc == n) goto finalize; // c is full...
        }
      }
    }
  }

#undef V

finalize:

  mad_free_tmp(U);
  mad_free_tmp(V);
  mad_free_tmp(S);

  // Return sorted indexes of columns to remove.
  return ivec_sort(c, nc, true);
}

int // Matrix preconditionning using SVD, return indexes of columns to remove.
mad_cmat_svdcnd(const cpx_t a[], idx_t c[], ssz_t m, ssz_t n,
               ssz_t N, num_t rcond, num_t s_[], num_t tol)
{
  assert(a && c);
  ssz_t mn = MIN(m,n);

  mad_alloc_tmp(cpx_t, U, m*m);
  mad_alloc_tmp(cpx_t, V, n*n);
  mad_alloc_tmp(num_t, S, mn );

  int info = mad_cmat_svd(a, U, S, V, m, n);
  if (info != 0) return -1;

  // Backup singular values.
  if (s_) mad_vec_copy(S, s_, mn);

  // N == 0 means to check for all singular values.
  if (N > mn || N <= 0) N = mn;

  // Tolerance on components similarity in V columns.
  tol = MAX(tol, DBL_EPSILON);

  // Tolerance on keeping singular values.
  rcond = MAX(rcond, DBL_EPSILON);

  // Number of columns to remove.
  idx_t nc = 0;

#define V(i,j) V[(i)*n+(j)]

  // Loop over increasing singular values.
  RFOR(i,mn,mn-N) {
    // Singular value is large, stop checking.
    if (S[i] > rcond*S[0]) break;

    // Loop over rows of V (i.e. columns of V^T)
    FOR(j,n-1) FOR(k,j+1,n) {
      num_t vj = cabs(V(j,i));

      // Proceed only significant component for this singular value.
      if (vj > 1e-4) {
        num_t vk  = cabs(V(k,i));
        num_t rat = fabs(vj-vk)/(vj+vk);

        // Discard column j with similar (or opposite) effect of column k > j.
        if (rat <= tol) {
          c[nc++] = j; // can hold duplicated indexes...
          if (nc == n) goto finalize; // c is full...
        }
      }
    }
  }

#undef V

finalize:

  mad_free_tmp(U);
  mad_free_tmp(V);
  mad_free_tmp(S);

  // Return sorted indexes of columns to remove.
  return ivec_sort(c, nc, true);
}

int // Matrix reconditionning using SVD.
mad_mat_pcacnd(const num_t a[], idx_t c[], ssz_t m, ssz_t n,
               ssz_t N, num_t rcond, num_t s_[])
{
  assert(a && c);
  ssz_t mn = MIN(m,n);

  mad_alloc_tmp(num_t, U, m*m);
  mad_alloc_tmp(num_t, V, n*n);
  mad_alloc_tmp(num_t, S, mn );
  mad_alloc_tmp(num_t, P, n  );

  int info = mad_mat_svd(a, U, S, V, m, n);
  if (info != 0) return -1;

  // Backup singular values.
  if (s_) mad_vec_copy(S, s_, mn);

  // N <= 0 means keep all columns.
  if (N > n || N <= 0) N = n;

  // Tolerance on keeping singular values.
  rcond = MAX(rcond, DBL_EPSILON);

  FOR(i,N) if (S[i] <= rcond*S[0]) { N=i; break; }

  // Compute projections on Principal Components, i.e. S V.
  mad_vec_abs(V, V, N*n);
  mad_mat_mul(S, V, P, 1, n, N);

  // Sort projections by ascending order.
  vec_sort(P, c, n);

  mad_free_tmp(U);
  mad_free_tmp(V);
  mad_free_tmp(S);
  mad_free_tmp(P);

  // Return sorted indexes of columns to remove.
  return ivec_sort(c, n-N, false);
}

int // Matrix reconditionning using SVD.
mad_cmat_pcacnd(const cpx_t a[], idx_t c[], ssz_t m, ssz_t n,
                ssz_t N, num_t rcond, num_t s_[])
{
  assert(a && c);
  ssz_t mn = MIN(m,n);

  mad_alloc_tmp(cpx_t, U, m*m);
  mad_alloc_tmp(cpx_t, V, n*n);
  mad_alloc_tmp(num_t, R, n*n);
  mad_alloc_tmp(num_t, S, mn );
  mad_alloc_tmp(num_t, P, n  );

  int info = mad_cmat_svd(a, U, S, V, m, n);
  if (info != 0) return -1;

  // Backup singular values.
  if (s_) mad_vec_copy(S, s_, mn);

  // N <= 0 means keep all columns.
  if (N > n || N <= 0) N = n;

  // Tolerance on keeping singular values.
  rcond = MAX(rcond, DBL_EPSILON);

  FOR(i,N) if (S[i] <= rcond*S[0]) { N=i; break; }

  // Compute projections on Principal Components, i.e. S V.
  mad_cvec_abs(V, R, N*n);
  mad_mat_mul (S, R, P, 1, n, N);

  // Sort projections by ascending order.
  vec_sort(P, c, n);

  mad_free_tmp(U);
  mad_free_tmp(V);
  mad_free_tmp(R);
  mad_free_tmp(S);
  mad_free_tmp(P);

  // Return sorted indexes of columns to remove.
  return ivec_sort(c, n-N, false);
}

static int // madx legacy code wrapper
madx_micado (const num_t a[], const num_t b[], num_t x[], ssz_t m, ssz_t n,
             ssz_t N, num_t tol, num_t r_[])
{
  /* copy buffers */
  mad_alloc_tmp(num_t, X   , n);
  mad_alloc_tmp(num_t, R   , m);
  /* working buffers */
  mad_alloc_tmp(idx_t, nx  , n);
  mad_alloc_tmp(idx_t, ny  , n);
  mad_alloc_tmp(num_t, ax  , m*n);
  mad_alloc_tmp(num_t, cinx, n);
  mad_alloc_tmp(num_t, xinx, m);
  mad_alloc_tmp(num_t, resx, m);
  mad_alloc_tmp(num_t, rho , 3*n);
  mad_alloc_tmp(num_t, ptop, n);
  mad_alloc_tmp(num_t, rmss, n);
  mad_alloc_tmp(num_t, xrms, n);
  mad_alloc_tmp(num_t, xptp, n);
  mad_alloc_tmp(num_t, xitr, n);

  mad_mat_trans(a, ax  , m, n);
  mad_vec_copy (b, xinx, m);
  mad_vec_fill (0, x   , n);

  int im=m, ic=n, iter=N, ifail=0;
  num_t rms=tol;

  micit_(X, R, nx, &rms, &im, &ic, &iter,
         /* working buffers */
         ny, ax, cinx, xinx, resx, rho, ptop, rmss, xrms, xptp, xitr, &ifail);

  // Re-order corrector strengths and save residues. Strengths are not minused!
  FOR(i,iter) x[i] = -X[nx[i]-1];
  if (r_) mad_vec_copy(R, r_, m);

  /* copy buffers */
  mad_free_tmp(X);
  mad_free_tmp(R);
  /* working buffers */
  mad_free_tmp(nx);
  mad_free_tmp(ny);
  mad_free_tmp(ax);
  mad_free_tmp(cinx);
  mad_free_tmp(xinx);
  mad_free_tmp(resx);
  mad_free_tmp(rho);
  mad_free_tmp(ptop);
  mad_free_tmp(rmss);
  mad_free_tmp(xrms);
  mad_free_tmp(xptp);
  mad_free_tmp(xitr);

  return iter;
}

/*
  Reference:
  B. Autin, and Y. Marti,
  "Closed Orbit Correction of A.G. Machines Using A Small Number of Magnets",
  CERN ISR-MA/73-17, 1973.
*/

int // Micado (from MAD9 + bug fixes)
mad_mat_nsolve(const num_t a[], const num_t b[], num_t x[], ssz_t m, ssz_t n,
               ssz_t N, num_t tol, num_t r_[])
{
  assert(a && b && x);

  // Micado notes: min_x || b - Ax ||_2 using N correctors amongst n
  // a: response matrix      [m x n]
  // b: vector of monitors   [m]
  // x: vector of correctors [n] (out)
  // r: residues             [m] (out)
  // N: number of correctors to use 0 < N <= n (out: actually used)

  mad_vec_fill(0, x, n);

  // No correctors.
  if (n == 0) return 0;

  if (tol < DBL_EPSILON) tol = DBL_EPSILON;

  // Checks if tolerance is already reached.
  { num_t e = mad_vec_norm(b, m) / m;
    if (e <= tol) return 0;
  }

  // N == 0 means to use all correctors
  if (N > n || N <= 0) N = n;

  // X-check with MAD-X Micado
  if (mad_use_madx_micado) {
    return madx_micado(a, b, x, m, n, N, tol, r_);
  }

  mad_alloc_tmp(num_t, A  , m*n);
  mad_alloc_tmp(num_t, B  , m);
  mad_alloc_tmp(num_t, X  , n);
  mad_alloc_tmp(num_t, R  , m);
  mad_alloc_tmp(num_t, sqr, n); // rho[k] = dot[k]^2/sqr[k]
  mad_alloc_tmp(num_t, dot, n);
  mad_alloc_tmp(idx_t, pvt, n);

  mad_vec_copy(a, A, m*n);
  mad_vec_copy(b, B, m);
  mad_vec_fill(0, X, n);
  mad_vec_fill(0, R, m);

#define A(i,j) A[(i)*n+(j)]

  // Box 3: Compute scalar products sqr[k] = A[k].A[k] and dot[k] = A[k].B.
  num_t sqrmin = 0;
  { num_t sum = 0;
    FOR(k,n) {
      pvt[k] = k;
      num_t hh = 0, gg = 0;
      FOR(i,m) {
        hh += A(i,k) * A(i,k);  // corrector effectiveness versus measured orbit
        gg += A(i,k) * B[i];    // corrector effectiveness versus target   orbit
      }
      sum   += hh;
      sqr[k] = hh;
      dot[k] = gg;
    }
    sqrmin = 1e-8 * sum / n;       // 1e-8 * average of correctors effectiveness
  }

  // Begin of iteration (l): loop over best-kick selection (i.e. A columns).
  FOR(k,N) {
    // Box 3: Search the columns not yet used for largest scaled change vector.
    { num_t maxChange = 0;
      idx_t changeIndex = -1;
      FOR(j,k,n) {
        if (sqr[j] > sqrmin) {        // criteria rho that minimize the residues
          num_t change = dot[j]*dot[j] / sqr[j];
          if (change > maxChange) {
            changeIndex = j;
            maxChange = change;
          }
        }
      }

      // Stop iterations if no suitable column are found.
      if (changeIndex < 0) { N=k; break; }

      // Move the column just found to next position.
      if (changeIndex > k) {
        num_t tr; idx_t ti;
        SWAP(sqr[k], sqr[changeIndex], tr);
        SWAP(dot[k], dot[changeIndex], tr);
        SWAP(pvt[k], pvt[changeIndex], ti);
        FOR(i,m) SWAP(A(i,k), A(i,changeIndex), tr);
      }
    }

    // Box 4: Compute beta, sigma, and vector u[k].
    num_t beta, hh;
    { hh = 0;
      FOR(i,k,m) hh += A(i,k) * A(i,k);
      num_t sigma = A(k,k) < 0 ? -sqrt(hh) : sqrt(hh);
      sqr[k] = -sigma; // saved for use in X[1..k] update
      A(k,k) += sigma;
      beta = 1 / (A(k,k) * sigma);
    }

    // Box 5: Transform remaining columns of A.
    FOR(j,k+1,n) {
      hh = 0;
      FOR(i,k,m) hh += A(i,k) * A(i,j);
      hh *= beta;
      FOR(i,k,m) A(i,j) -= A(i,k) * hh;
    }

    // Box 6: Transform vector b.
    hh = 0;
    FOR(i,k,m) hh += A(i,k) * B[i];
    hh *= beta;
    FOR(i,k,m) B[i] -= A(i,k) * hh;

    // Box 3: Update scalar products sqr[j]=A[j]*A[j] and dot[j]=A[j]*b.
    FOR(j,k+1,n) {
      sqr[j] -= A(k,j) * A(k,j);
      dot[j] -= A(k,j) * B[k];
    }

    // Box 7: Recalculate solution vector x. Here, sqr[1..k] = -sigma[1..k].
    RFOR(i,k+1) {
      X[i] = B[i];
      FOR(j,i+1,k+1) X[i] -= A(i,j) * X[j];
      X[i] /= sqr[i];
    }

    // Box 8: Compute original residual vector by backward transformation.
    mad_vec_copy(B, R, m);
    RFOR(j,k+1) {
      R[j] = hh = 0;
      FOR(i,j,m) hh += A(i,j) * R[i];
      hh /= sqr[j] * A(j,j);
      FOR(i,j,m) R[i] += A(i,j) * hh;
    }

    // Box 9: Check for convergence.
    num_t e = mad_vec_norm(R, m) / m;
    if (e <= tol) { N=k+1; break; }
  }

#undef A

  // Re-order corrector strengths and save residues.
  FOR(i,N) x[pvt[i]] = X[i];
  if (r_) mad_vec_copy(R, r_, m);

  mad_free_tmp(A);
  mad_free_tmp(B);
  mad_free_tmp(X);
  mad_free_tmp(R);
  mad_free_tmp(sqr);
  mad_free_tmp(dot);
  mad_free_tmp(pvt);

  return N;
}

// -- Survey Misalignments ----------------------------------------------------o

#define N 3

// Fast computation of Rbar and Tbar for misalignments at element exit (forward)
// Tbar = W^t (RV+T-V)
// Rbar = W^t R W

#include "mad_cst.h"

void
mad_mat_rtbar (num_t Rb[NN],       num_t Tb[N], num_t el, num_t ang, num_t tlt,
         const num_t R_[NN], const num_t T [N])
{
  assert(Rb && Tb && T);

  if (fabs(ang) < mad_cst_MINANG) {           // -- straight ------------------o
    if (R_) {
      num_t Ve[N] = {0, 0, el};               // We = I
      mad_mat_mul (R_, Ve, Tb, N, 1, N);
      mad_vec_sub (Tb, Ve, Tb, N);
      mad_vec_add (Tb, T , Tb, N);            // Tb = R*Ve + T - Ve
      mad_vec_copy(R_,     Rb, NN);           // Rb = R
    } else { // R = I
      mad_vec_copy(T, Tb, N);                 // Tb = T
      mad_mat_eye (Rb, 1, N, N, N);           // Rb = I
    }

  } else {                                    // -- curved --------------------o
    num_t rho = el/ang;
    num_t Ve[N] = {rho*(cos(ang)-1), 0, rho*sin(ang)};
    num_t We[NN]; mad_mat_roty(We, -ang);
    num_t Wt[NN];

    // tilt
    if (fabs(tlt) >= mad_cst_MINANG) {
      num_t Rz[NN]; mad_mat_rotz(Rz, tlt);
      mad_mat_mul (Rz, Ve, Ve, N, 1, N);      // Ve = Rz*Ve
      mad_mat_mul (Rz, We, Wt, N, N, N);
      mad_mat_mult(Wt, Rz, We, N, N, N);      // We = Rz*We*Rz:t()
    }

    if (R_) {
      num_t Vt[N];
      mad_mat_mul (R_, Ve, Vt, N, 1, N);
      mad_vec_sub (Vt, Ve, Vt, N);
      mad_vec_add (Vt, T , Vt, N);
      mad_mat_tmul(We, Vt, Tb, N, 1, N);      // Tb = We:t()*(R*Ve + T - Ve)
      mad_mat_tmul(We, R_, Wt, N, N, N);
      mad_mat_mul (Wt, We, Rb, N, N, N);      // Rb = We:t()*R*We
    } else { // R = I
      mad_mat_tmul(We, T , Tb, N, 1, N);      // Tb = We:t()*T
      mad_mat_eye (Rb,      1, N, N, N);      // Rb = I
    }
  }
}

#undef N
