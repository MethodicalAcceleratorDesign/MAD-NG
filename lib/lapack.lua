local oss  = jit.os
local ffi  = require 'ffi'
local lapack = ffi.load("lib/lapack/liblapack-" .. oss .. ".so")

--[[
LAPACK is the default method for solving dense numerical matrices. When the
matrix is square and nonsingular the routines dgesv, dlange, and dgecon are used
for real matrices and zgesv, zlange, and zgecon for complex matrices. When the
matrix is nonsquare or singular dgelss is used for real matrices and zgelss for
complex matrices. More information on LAPACK is available in the references.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. For generalized eigenvalues the routine dggev is used for real
matrices and zggev for complex matrices. When the matrix is unsymmetric, dgeev
is used for real matrices and zgeev is used for complex matrices. For symmetric
matrices, dsyevr is used for real matrices and zheevr is used for complex matrices.
]]

--[[
Type convertion:
Fortran             C
CHARACTER           char*
INTEGER             int*
DOUBLE PRECISION    double*

Transpose:
C       -> Row Major
Fortran -> Col Major (transpose of C)
'T' (real) or 'C' (complex conjugate)
]]

ffi.cdef[[
// -----
// C <- alpha * A * B + beta * C with A[m x p], B[p x n] and C[m x n]                                       
// if beta ~= 0.0, C must be initialized
// slow functions, replaced by faster C matmul
// -----
void dgemm_(const char *tA, const char *tB, const int *m, const int *n, const int *p,
            const double  *alpha, const double  *A, const int *lda,
                                  const double  *B, const int *ldb,
            const double  *beta ,       double  *C, const int *ldc);
void zgemm_(const char *tA, const char *tB, const int *m, const int *n, const int *p,
            const complex *alpha, const complex *A, const int *lda,
                                  const complex *B, const int *ldb,
            const complex *beta ,       complex *C, const int *ldc);

// -----
// inverse A[m x n]
// A <- A^-1, ipiv[min(m,n)] <- pivot indices
// info < 0 -> invalid input, info > 0 -> singular matrix
// sizes: lwork=-1, work[0]=optimal size
// -----
void dgetrf_(const int* m, const int *n, double* A, const int* lda, int* ipiv, int* info);
void dgetri_(const int* n, double* A, const int* lda, int* ipiv, double* work, int* lwork, int* info);

void zgetrf_(const int* m, const int *n, complex* A, const int* lda, int* ipiv, int* info);
void zgetri_(const int* n, complex* A, const int* lda, int* ipiv, complex* work, int* lwork, int* info);

// -----
// solve A * X = B with A[n x n] and B[n x nrhs]
// A <- L * U, ipiv[n] <- pivot indices, B <- X (if info = 0)
// info < 0 -> invalid input, info > 0 -> singular matrix
// -----
void dgesv_(const int *n, const int *nrhs, double  *A, const int *lda,
                                int *IPIV, double  *B, const int *ldb, int *info);
void zgesv_(const int *n, const int *nrhs, complex *A, const int *lda,
                                int *IPIV, complex *B, const int *ldb, int *info);

// -----
// solve A * X = B with A[n x n], B[n x nrhs] and X[n x nrhs] to machine precision
// A <- L * U, ipiv[n] <- pivot indices, X <- solution
// info < 0 -> invalid input, info > 0 -> singular matrix
// sizes: AF[n x n], R[n], C[n], ferr[nrhs], berr[nrhs], work[4n], iwork[n]
// -----
void dgesvx_(const char *fact, const char *trans, const int *n, const int *nrhs,
             double *A, const int *lda, double *AF, const int *ldaf, int *ipiv,
             const char *equed, double *R, double *C,
             double *B, const int *ldb,
             double *X, const int *ldx, double *rcond, double *ferr, double *berr,
             double *work, int *iwork, int *info);
void zgesvx_(const char *fact, const char *trans, const int *n, const int *nrhs,
             complex *A, const int *lda, complex *AF, const int *ldaf, int *ipiv,
             const char *equed, double *R, double *C,
             complex *B, const int *ldb,
             complex *X, const int *ldx, double *rcond, double *ferr, double *berr,
             complex *work, double *rwork, int *info);

// -----
// Solve A * X = B with A[m x n], B[m x nrhs] and X[m x nrhs]: search min | b - Ax | using SVD
// A <- right singular vectors, B <- X (if info = 0), S <- singular value in decreasing order
// rcond: S(i) <= rcond*S(1) are retained. rcond < 0 -> use machine precision
// rank <- rank of A (i.e. number of SV)
// info < 0 -> invalid input, info > 0 -> singular matrix
// R sizes: work[lwork], lwork >= 3*min(m,n) + max( 2*min(m,n), max(m,n), nrhs )
// C sizes: work[lwork], lwork >= 2*min(m,n) + max(m,n,nrhs), rwork[5*min(m,n)]
// -----
void dgelss_(const int *m, const int *n, const int *nrhs,
             double *A, const int *lda, double *B, const int *ldb,
             double *S, const double *rcond, int *rank,
             double *work, const int *lwork, int *info);
void zgelss_(const int *m, const int *n, const int *nrhs,
             complex *A, const int *lda, complex *B, const int *ldb,
             double *S, const double *rcond, int *rank,
             complex *work, const int *lwork, double *rwork, int *info);

// Eigen values/vectors
// void dggev_
// void zggev_
]]

return lapack