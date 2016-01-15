/*
 o----------------------------------------------------------------------------o
 |
 | Matrix module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

#include <math.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"

// --- implementation --------------------------------------------------------o

#define CHKR     assert( r )
#define CHKXY    assert( x && y )
#define CHKXR    assert( x && r )
#define CHKYR    assert( y && r )
#define CHKXYR   assert( x && y && r )
#define CHKXRX   assert( x && r && x != r)
#define CHKYRY   assert( y && r && y != r)
#define CHKXYRX  assert( x && y && r && x != r )
#define CHKXYRY  assert( x && y && r && y != r )
#define CHKXYRXY assert( x && y && r && x != r && y != r )

#define CNUM(a) cnum_t a = (* (cnum_t*) & (num_t[2]) { (a##_re), (a##_im) })

// [m x n] = [m x p] * [p x n]
// naive implementation (not vectorized)
/* #define MMUL \
  for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
*/

// [m x n] = [m x p] * [p x n]
// portable vectorized matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
// get xN speed-up factor compared to dgemm from openblas and lapack...
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

#define MUL() { \
  if (m == 1 && n == 1) VMULV \
  else      if (m == 1) VMUL \
  else      if (n == 1) MULV \
  else      if (p == 1) KMUL \
  else                  MMUL \
} \

// -----

// r = trace([p x m]^t * [p x n])
// Frobenius inner product
#define DOT(C) { \
  *r = 0; \
  if (m == 1 && n == 1) { \
    for (size_t i=0; i < p; i++) \
      *r += C(x[i]) * y[i]; \
  } else { \
    size_t mn = MIN(m,n); \
    for (size_t i=0; i < mn; i++) \
      for (size_t k=0; k < p; k++) \
        *r += C(x[k*m+i]) * y[k*n+i]; \
  } \
} \

// [m x n] transpose, in place only for square matrix
// inplace for non-square matrix could use FFTW...
// see http://www.fftw.org/faq/section3.html#transpose
#define TRANS(C) { \
  if (m == 1 || n == 1) { \
    size_t mn = m*n; \
    for (size_t i=0; i < mn; i++) \
      r[i] = C(x[i]); \
  } else if ((const void*)x == (const void*)r) { \
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

#define CPY(OP) { \
  for (size_t i=0; i<m; i++) \
  for (size_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x[i*ldx+j]; \
} \

#define SET(OP) { \
  for (size_t i=0; i<m; i++) \
  for (size_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x; \
} \

#define DIAG(OP) { \
  for (size_t i=0; i<MIN(m,n); i++) \
    r[i*ldr+i] OP##= x; \
} \

// --- mat

void mad_mat_ident(num_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; num_t x = 0; SET(); x = 1; DIAG(); }

void mad_mat_set(num_t x, num_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; SET(); }

void mad_mat_cpy(const num_t x[], num_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXRX; CPY(); }

void mad_mat_cpym(const num_t x[], cnum_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXR; CPY(); }

void mad_mat_trans (const num_t x[], num_t r[], size_t m, size_t n)
{ CHKXR; TRANS(); }

num_t mad_mat_dot (const num_t x[], const num_t y[], size_t m, size_t n, size_t p)
{ CHKXY; num_t r_, *r=&r_; DOT(); return *r; }

void mad_mat_dotm (const num_t x[], const cnum_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(); }

void mad_mat_mul (const num_t x[], const num_t y[], num_t r[], size_t m, size_t n, size_t p)
{ CHKXYRXY; MUL(); }

void mad_mat_mulm (const num_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRY; MUL(); }

// -- cmat

void mad_cmat_ident(cnum_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; cnum_t x = 0; SET(); x = 1; DIAG(); }

void mad_cmat_set(num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; CNUM(x); SET(); }

void mad_cmat_cpy(const cnum_t x[], cnum_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXRX; CPY(); }

void mad_cmat_trans (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{ CHKXR; TRANS(); }

void mad_cmat_ctrans (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{ CHKXR; TRANS(conj); }

void mad_cmat_dot (const cnum_t x[], const cnum_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(conj); }

void mad_cmat_dotm (const cnum_t x[], const num_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXY; DOT(conj); }

void mad_cmat_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRXY; MUL(); }

void mad_cmat_mulm (const cnum_t x[], const num_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRX; MUL(); }

// -- lapack -----------------------------------------------------------------o

/*
LAPACK is the default method for solving dense numerical matrices. When the
matrix is square and nonsingular the routines dgesv and zgesv are used
otherwise routines dgelsy and zgelsy are used.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. For generalized eigenvalues the routine dggev zggev are used.
See also [d,z]geev and [d,z]syevr, zheevr.
*/

// -----
// solve A * X = B with A[n x n], B[n x nrhs] and X[n x nrhs]: search min | b - Ax | using LU
// -----
void dgesv_ (const int *n, const int *nrhs, num_t  *A, const int *lda,
                                 int *IPIV, num_t  *B, const int *ldb, int *info);
void zgesv_ (const int *n, const int *nrhs, cnum_t *A, const int *lda,
                                 int *IPIV, cnum_t *B, const int *ldb, int *info);

// -----
// Solve A * X = B with A[m x n], B[m x nrhs] and X[m x nrhs]: search min | b - Ax | using QR
// -----
void dgelsy_ (const int *m, const int *n, const int *nrhs,
              num_t *A, const int *lda, num_t *B, const int *ldb,
              int *jpvt, const num_t *rcond, int *rank,
              num_t *work, const int *lwork, int *info);
void zgelsy_ (const int *m, const int *n, const int *nrhs,
              cnum_t *A, const int *lda, cnum_t *B, const int *ldb,
              int *jpvt, const num_t *rcond, int *rank,
              cnum_t *work, const int *lwork, num_t *rwork, int *info);

// Eigen values/vectors
// void dggev_
// void zggev_

// -- inverse ----------------------------------------------------------------o

int
mad_mat_invn (const num_t y[], num_t x, num_t r[], size_t m, size_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  alloc_tmp(num_t, u, n*n);
  mad_mat_ident(u, n, n, n);
  int rank = mad_mat_div(u, y, r, n, m, n, rcond);
  free_tmp(u);
  if (x != 1.0) mad_vec_muln(r, x, r, m*n);
  return rank;
}

int
mad_mat_invc (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  alloc_tmp(num_t, u, n*n);
  mad_mat_ident(u, n, n, n);
  alloc_tmp(num_t, t, m*n);
  int rank = mad_mat_div(u, y, t, n, m, n, rcond);
  free_tmp(u);
  if (x_re != 1.0 || x_im != 0.0) mad_vec_mulc(t, x_re, x_im, r, m*n);
  free_tmp(t);
  return rank;
}

int
mad_cmat_invn (const cnum_t y[], num_t x, cnum_t r[], size_t m, size_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  alloc_tmp(cnum_t, u, n*n);
  mad_cmat_ident(u, n, n, n);
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
  free_tmp(u);
  if (x != 1.0) mad_cvec_muln(r, x, r, m*n);
  return rank;
}

int
mad_cmat_invc (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, num_t rcond)
{
  CHKYR; // compute U:[n x n]/Y:[m x n]
  alloc_tmp(cnum_t, u, n*n);
  mad_cmat_ident(u, n, n, n);
  int rank = mad_cmat_div(u, y, r, n, m, n, rcond);
  free_tmp(u);
  if (x_re != 1.0 || x_im != 0.0)  mad_cvec_mulc(r, x_re, x_im, r, m*n);
  return rank;
}

// -- divide -----------------------------------------------------------------o

// note:
// X/Y => X * Y^-1 => [m x p] * [p x n] => X:[m x p], Y:[n x p]
// X/Y => X * Y^-1 => (Y'^-1 * X')' => A=Y' and B=X'
// Solving A*X=B => X = A^-1 B = (B'/A')' (col-major!)
//    with Y':[p x n] = A:[M=p x N=n],
//    and  X':[p x m] = B:[M=p x NRHS=m], ipiv:[N]

int
mad_mat_div (const num_t x[], const num_t y[], num_t r[], size_t m, size_t n, size_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nn=n, nm=m, np=p;
  alloc_tmp(num_t, a, n*p);
  mad_vec_cpy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) { 
    int ipiv[n];
    mad_vec_cpy(x, r, m*p);
    dgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  num_t sz;
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  alloc_tmp(num_t, rr, ldb*nm);
  mad_mat_cpy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  dgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, &info); // query
  alloc_tmp(num_t, wk, lwork=sz);
  dgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, &info); // compute
  mad_mat_cpy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  free_tmp(wk); free_tmp(rr); free_tmp(a);

  if (info < 0) fatal("invalid input argument");
  if (info > 0) fatal("unexpect lapack error");

  return rank;
}

int
mad_mat_divm (const num_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nn=n, nm=m, np=p;
  alloc_tmp(cnum_t, a, n*p);
  mad_cvec_cpy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) { 
    int ipiv[n];
    mad_vec_cpyv(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  alloc_tmp(cnum_t, rr, ldb*nm);
  mad_mat_cpym(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_cpy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  free_tmp(wk); free_tmp(rr); free_tmp(a);

  if (info < 0) fatal("invalid input argument");
  if (info > 0) fatal("unexpect lapack error");

  return rank;
}

int
mad_cmat_div (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nn=n, nm=m, np=p;
  alloc_tmp(cnum_t, a, n*p);
  mad_cvec_cpy(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) { 
    int ipiv[n];
    mad_cvec_cpy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  alloc_tmp(cnum_t, rr, ldb*nm);
  mad_cmat_cpy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_cpy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  free_tmp(wk); free_tmp(rr); free_tmp(a);

  if (info < 0) fatal("invalid input argument");
  if (info > 0) fatal("unexpect lapack error");

  return rank;
}

int
mad_cmat_divm (const cnum_t x[], const num_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
{
  CHKXYR;
  int info=0;
  const int nn=n, nm=m, np=p;
  alloc_tmp(cnum_t, a, n*p);
  mad_vec_cpyv(y, a, n*p);

  // square system (y is square, n == p), use LU decomposition
  if (n == p) { 
    int ipiv[n];
    mad_cvec_cpy(x, r, m*p);
    zgesv_(&np, &nm, a, &np, ipiv, r, &np, &info);
    if (!info) return free_tmp(a), n;
  }

  // non-square system or singular square system, use QR or LQ factorization
  cnum_t sz;
  num_t rwk[2*nn];
  int rank, ldb=MAX(nn,np), lwork=-1; // query for optimal size
  int JPVT[nn]; memset(JPVT, 0, sizeof JPVT);
  alloc_tmp(cnum_t, rr, ldb*nm);
  mad_cmat_cpy(x, rr, m, p, p, ldb); // input strided copy [M x NRHS]
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank, &sz, &lwork, rwk, &info); // query
  alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zgelsy_(&np, &nn, &nm, a, &np, rr, &ldb, JPVT, &rcond, &rank,  wk, &lwork, rwk, &info); // compute
  mad_cmat_cpy(rr, r, m, n, ldb, n); // output strided copy [N x NRHS]
  free_tmp(wk); free_tmp(rr); free_tmp(a);

  if (info < 0) fatal("invalid input argument");
  if (info > 0) fatal("unexpect lapack error");

  return rank;
}

