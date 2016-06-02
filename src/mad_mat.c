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

#define NO(a) 
#define ID(a) a

#define CNUM(a) cnum_t a = (* (cnum_t*) & (num_t[2]) { MKNAME(a,_re), MKNAME(a,_im) })

// [m x n] = [m x p] * [p x n]
// naive implementation (not vectorized)
/* #define MMUL \
  for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    }
*/

// [m x n] = [m x p] * [p x n]
// portable vectorized gerenal matrix-matrix multiplication
// loop unroll + vectorized on SSE2 (x2), AVX & AVX2 (x4), AVX-512 (x8)
// get xN speed-up factor compared to dgemm from openblas and lapack...
#define MMUL { /* mat * mat */ \
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
}

// n==1: [m x 1] = [m x p] * [p x 1]
#define MULV /* mat * vec */ \
  for (size_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (size_t k=0; k < p; k++) \
      r[i] += x[i*p+k] * y[k]; \
  }

// m==1: [1 x n] = [1 x p] * [p x n]
#define VMUL /* vec * mat */ \
  for (size_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (size_t k=0; k < p; k++) \
      r[j] += x[k] * y[k*n+j]; \
  }

// p==1: [m x n] = [m x 1] * [1 x n]
#define KMUL /* outer */ \
  for (size_t i=0; i < m; i++) { \
  for (size_t j=0; j < n; j++) \
    r[i*n+j] = x[i] * y[j]; \
  }

// m==1, n==1: [1 x 1] = [1 x p] * [p x 1]
#define VMULV /* inner */ \
  { r[0] = 0; \
    for (size_t k=0; k < p; k++) \
      r[0] += x[k] * y[k]; \
  }

#define MUL() { \
  if (m == 1 && n == 1) VMULV \
  else      if (m == 1) VMUL \
  else      if (n == 1) MULV \
  else      if (p == 1) KMUL \
  else                  MMUL \
}

// -----

// r = trace([p x m]^t * [p x n])
// Frobenius inner product
#define DOT(C) { \
  if (m == 1 && n == 1) { \
    for (size_t i=0; i < p; i++) \
      *r += C(x[i]) * y[i]; \
  } else { \
    size_t mn = MIN(m,n); \
    for (size_t i=0; i < mn; i++) \
      for (size_t k=0; k < p; k++) \
        *r += C(x[k*m+i]) * y[k*n+i]; \
  } \
}

// [m x n] transpose
#define TRANS(C,T) { \
  if (m == 1 || n == 1) { \
    size_t mn = m*n; \
    for (size_t i=0; i < mn; i++) \
      r[i] = C(x[i]); \
  } else if ((const void*)x != (const void*)r) { \
    for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) \
      r[j*m+i] = C(x[i*n+j]); \
  } else if (m == n) { \
    for (size_t i=0; i < m; i++) \
    for (size_t j=i; j < n; j++) { \
      T t = C(r[j*m+i]); \
      r[j*m+i] = C(r[i*n+j]); \
      r[i*n+j] = t; \
    } \
  } else { \
    mad_alloc_tmp(T, t, m*n); \
    for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) \
      t[j*m+i] = C(x[i*n+j]); \
    memcpy(r, t, m*n*sizeof(T)); \
    mad_free_tmp(t); \
  } \
}

#define CPY(OP) { \
  for (size_t i=0; i<m; i++) \
  for (size_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x[i*ldx+j]; \
}

#define SET(OP) { \
  for (size_t i=0; i<m; i++) \
  for (size_t j=0; j<n; j++) \
    r[i*ldr+j] OP##= x; \
}

#define DIAG(OP) { \
  for (size_t i=0; i<MIN(m,n); i++) \
    r[i*ldr+i] OP##= x; \
}

// --- mat

void mad_mat_ident(num_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; num_t x = 0; SET(); x = 1; DIAG(); }

void mad_mat_fill(num_t x, num_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; SET(); }

void mad_mat_copy(const num_t x[], num_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXRX; CPY(); }

void mad_mat_copym(const num_t x[], cnum_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXR; CPY(); }

void mad_mat_trans (const num_t x[], num_t r[], size_t m, size_t n)
{ CHKXR; TRANS(,num_t); }

num_t mad_mat_dot (const num_t x[], const num_t y[], size_t m, size_t n, size_t p)
{ CHKXY; num_t r_=0, *r=&r_; DOT(ID); return *r; }

cnum_t mad_mat_dotm (const num_t x[], const cnum_t y[], size_t m, size_t n, size_t p)
{ CHKXY; cnum_t r_=0, *r = &r_; DOT(ID); return *r; }

void mad_mat_dotm_r (const num_t x[], const cnum_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; *r=0; DOT(ID); }

void mad_mat_mul (const num_t x[], const num_t y[], num_t r[], size_t m, size_t n, size_t p)
{ CHKXYRXY; MUL(); }

void mad_mat_mulm (const num_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRY; MUL(); }

void mad_mat_center (const num_t x[], num_t r[], size_t m, size_t n)
{
  mad_alloc_tmp(num_t, u, m);
  for (size_t i=0; i < m; i++) { u[i]  = 0;
  for (size_t j=0; j < n; j++)   u[i] += x[i*n+j]; }
  for (size_t i=0; i < m; i++)
  for (size_t j=0; j < n; j++) r[i*n+j] = x[i*n+j] - u[i];
  mad_free_tmp(u);
}

// -- cmat

void mad_cmat_ident(cnum_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; cnum_t x = 0; SET(); x = 1; DIAG(); }

void mad_cmat_fill(cnum_t x, cnum_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; SET(); }

void mad_cmat_fill_r(num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, size_t ldr)
{ CHKR; CNUM(x); SET(); }

void mad_cmat_copy(const cnum_t x[], cnum_t r[], size_t m, size_t n, size_t ldx, size_t ldr)
{ CHKXRX; CPY(); }

void mad_cmat_trans (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{ CHKXR; TRANS(,cnum_t); }

void mad_cmat_ctrans (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{ CHKXR; TRANS(conj,cnum_t); }

cnum_t mad_cmat_dot (const cnum_t x[], const cnum_t y[], size_t m, size_t n, size_t p)
{ CHKXY; cnum_t r_=0, *r = &r_; DOT(conj); return *r; }

void mad_cmat_dot_r (const cnum_t x[], const cnum_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; *r=0; DOT(conj); }

cnum_t mad_cmat_dotm (const cnum_t x[], const num_t y[], size_t m, size_t n, size_t p)
{ CHKXY; cnum_t r_=0, *r = &r_; DOT(conj); return *r; }

void mad_cmat_dotm_r (const cnum_t x[], const num_t y[], cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; *r=0; DOT(conj); }

void mad_cmat_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRXY; MUL(); }

void mad_cmat_mulm (const cnum_t x[], const num_t y[], cnum_t r[], size_t m, size_t n, size_t p)
{ CHKXYRX; MUL(); }

void mad_cmat_center (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{
  mad_alloc_tmp(cnum_t, u, m);
  for (size_t i=0; i < m; i++) { u[i]  = 0;
  for (size_t j=0; j < n; j++)   u[i] += x[i*n+j]; }
  for (size_t i=0; i < m; i++)
  for (size_t j=0; j < n; j++) r[i*n+j] = x[i*n+j] - u[i];
  mad_free_tmp(u);
}

// -- lapack -----------------------------------------------------------------o

/*
LAPACK is the default method for solving dense numerical matrices. When the
matrix is square and nonsingular the routines dgesv and zgesv are used
otherwise routines dgelsy and zgelsy are used.

LAPACK is the default method for computing the entire set of singluar values
and singular vectors. For generalized SVD the routines dgesdd and zgesdd are used.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. For generalized eigenvalues the routines dggev and zggev are used.
*/

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
void dggev_ (str_t jobvl, str_t jobvr, const int *n, num_t A[], const int *lda,
             num_t WR[], num_t WI[],
             num_t VL[], const int *ldvl, num_t VR[], const int *ldvr,
             num_t work[], int *lwork, int *info);
void zggev_ (str_t jobvl, str_t jobvr, const int *n, cnum_t A[], const int *lda,
             cnum_t W[], cnum_t VL[], const int *ldvl, cnum_t VR[], const int *ldvr,
             cnum_t work[], int *lwork, num_t rwork[], int *info);

// -- inverse ----------------------------------------------------------------o

int
mad_mat_invn (const num_t y[], num_t x, num_t r[], size_t m, size_t n, num_t rcond)
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
mad_mat_invc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, num_t rcond)
{ CNUM(x); return mad_mat_invc(y, x, r, m, n, rcond); }

int
mad_mat_invc (const num_t y[], cnum_t x, cnum_t r[], size_t m, size_t n, num_t rcond)
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
mad_cmat_invn (const cnum_t y[], num_t x, cnum_t r[], size_t m, size_t n, num_t rcond)
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
mad_cmat_invc (const cnum_t y[], cnum_t x, cnum_t r[], size_t m, size_t n, num_t rcond)
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
mad_cmat_invc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t m, size_t n, num_t rcond)
{ CNUM(x); return mad_cmat_invc(y, x, r, m, n, rcond); }

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
mad_mat_divm (const num_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
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
mad_cmat_div (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
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
mad_cmat_divm (const cnum_t x[], const num_t y[], cnum_t r[], size_t m, size_t n, size_t p, num_t rcond)
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
mad_mat_svd (const num_t x[], num_t u[], num_t s[], num_t v[], size_t m, size_t n)
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
mad_cmat_svd (const cnum_t x[], cnum_t u[], num_t s[], cnum_t v[], size_t m, size_t n)
{
  assert( x && u && s && v );
  int info=0;
  const int nm=m, nn=n;

  cnum_t sz;
  int lwork=-1;
  int iwk[8*MIN(m,n)];
  size_t rwk_sz = MIN(m,n) * MAX(5*MIN(m,n)+7, 2*MAX(m,n)+2*MIN(m,n)+1);
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
mad_mat_eigen (const num_t x[], cnum_t w[], num_t vl[], num_t vr[], size_t n)
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
  dggev_("V", "V", &nn, ra, &nn, wr, wi, vl, &nn, vr, &nn, &sz, &lwork, &info); // query
  mad_alloc_tmp(num_t, wk, lwork=sz);
  dggev_("V", "V", &nn, ra, &nn, wr, wi, vl, &nn, vr, &nn,  wk, &lwork, &info); // compute
  mad_vec_cvec(wi, wr, w, n);
  mad_free_tmp(wk); mad_free_tmp(ra);
  mad_free_tmp(wi); mad_free_tmp(wr); 
  mad_mat_trans(vl, vl, n, n);
  mad_mat_trans(vr, vr, n, n);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("eigen failed to compute all eigenvalues");

  return info;
}

int
mad_cmat_eigen (const cnum_t x[], cnum_t w[], cnum_t vl[], cnum_t vr[], size_t n)
{
  assert( x && w && vl && vr );
  int info=0;
  const int nn=n;

  cnum_t sz;
  int lwork=-1;
  mad_alloc_tmp(num_t, rwk, 2*n);
  mad_alloc_tmp(cnum_t, ra, n*n);
  mad_cmat_trans(x, ra, n, n);
  zggev_("V", "V", &nn, ra, &nn, w, vl, &nn, vr, &nn, &sz, &lwork, rwk, &info); // query
  mad_alloc_tmp(cnum_t, wk, lwork=creal(sz));
  zggev_("V", "V", &nn, ra, &nn, w, vl, &nn, vr, &nn,  wk, &lwork, rwk, &info); // compute
  mad_free_tmp(wk); mad_free_tmp(ra); mad_free_tmp(rwk);
  mad_cmat_trans(vl, vl, n, n);
  mad_cmat_trans(vr, vr, n, n);

  if (info < 0) error("invalid input argument");
  if (info > 0) warn ("eigen failed to compute all eigenvalues");

  return info;
}

// -- FFT ---------------------------------------------------------------------o

#include <fftw3.h>

void // x [m x n] -> r [m, n/2+1]
mad_mat_fft (const num_t x[], cnum_t r[], size_t m, size_t n)
{
  CHKXR;
  mad_alloc_tmp(cnum_t, cx, m*n);
  mad_vec_copyv(x, cx, m*n);
  mad_cmat_fft(cx, r, m, n);
  mad_free_tmp(cx);
}

void // x [m x n] -> r [m, n/2+1]
mad_mat_rfft (const num_t x[], cnum_t r[], size_t m, size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_r2c_2d(m, n, (num_t*)x, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void // x [m x n/2+1] -> r [m x n]
mad_mat_irfft (const cnum_t x[], num_t r[], size_t m, size_t n)
{
  CHKXR;
  size_t nn = m*(n/2+1);
  mad_alloc_tmp(cnum_t, cx, nn);
  mad_cvec_copy(x, cx, nn);
  fftw_plan p = fftw_plan_dft_c2r_2d(m, n, cx, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_free_tmp(cx);
  mad_vec_muln(r, 1.0/(m*n), r, m*n);
}

void
mad_cmat_fft (const cnum_t x[], cnum_t r[], size_t m, size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_2d(m, n, (cnum_t*)x, r, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cmat_ifft(const cnum_t x[], cnum_t r[], size_t m, size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_2d(m, n, (cnum_t*)x, r, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_cvec_muln(r, 1.0/(m*n), r, m*n);
}

// -- NFFT --------------------------------------------------------------------o
