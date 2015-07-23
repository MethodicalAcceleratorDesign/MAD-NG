local oss  = jit.os
local ffi  = require 'ffi'
local lapack = ffi.load("lib/lapack/liblapack-" .. oss .. ".so")

--[[
LAPACK is the default method for solving dense numerical matrices. When the
matrix is square and nonsingular the routines dgesv, dlange, and dgecon are used
for real matrices and zgesv, zlange, and zgecon for complex matrices. When the
matrix is nonsquare or singular dgelss is used for real matrices and zgelss for
complex matrices. More information on LAPACK is available in the references.

For dense matrices the Cholesky method uses LAPACK functions such as dpotrf and
dpotrs for real matrices and zpotrf and zpotrs for complex matrices. For sparse
matrices the Cholesky method uses the "TAUCS" library.

For eigenvalue computation when the input is an × matrix of machine numbers and
the number of eigenvalues requested is less than 20% of  an Arnoldi method is
used. Otherwise, if the input matrix is numerical then a method using LAPACK
routines is used.

LAPACK is the default method for computing the entire set of eigenvalues and
eigenvectors. When the matrix is unsymmetric, dgeev is used for real matrices
and zgeev is used for complex matrices. For symmetric matrices, dsyevr is used
for real matrices and zheevr is used for complex matrices. For generalized
eigenvalues the routine dggev is used for real matrices and zggev for complex
matrices.

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
"T" (real) or "C" (complex)

Example: C = alpha A * B + beta C

SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
.. Scalar Arguments ..
DOUBLE PRECISION ALPHA,BETA
INTEGER K,LDA,LDB,LDC,M,N
CHARACTER TRANSA,TRANSB
.. Array Arguments ..
DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

void dgemm_(const char* TRANSA, const char* TRANSB, const int* M,
            const int* N, const int* K, const double* alpha, const double* A,
            const int* LDA, const double* B, const int* LDB, const double* beta,
            double* C, const int* LDC);

Use case: degmm_('T', 'T', M, N, K, 1.0, A, K, B, N, 0.0, C, N)
]]

ffi.cdef[[
// -----
// C <- alpha * A * B + beta * C with A[m x p], B[p x n] and C[m x n]                                       
// if beta ~= 0.0, C must be initialized
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
// solve A * X = B with A[n x n] and B[n x nrhs]
// A <- L * U, ipiv[n] <- pivot indices, B <- X (if info = 0)
// info < 0 -> invalid input, info > 0 -> singular matrix
// -----
void dgesv_(const int *n, const int *nrhs, double  *A, const int *lda,
                                int *ipiv, double  *B, const int *ldb, int *info);
void zgesv_(const int *n, const int *nrhs, complex *A, const int *lda,
                                int *ipiv, complex *B, const int *ldb, int *info);

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

/*
          2: DGELSY 
          3: DGELSS 
          4: DGELSD 
*/

// M x N
void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
            double *A, const int *lda, double *b, const int *ldb, double *work,
            const int *lwork, int *info);
//void zgels_

// Driver M x N: min | b - Ax | using SVD
//void dgelss_
//void zgelss_


// SVD + driver
void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N,
             double* A, const int* lda, double* S, double* U, const int* ldu,
             double* VT, const int* ldvt, double* work,const int* lwork,
             const int* info);
//void zgesvd_
//void dgesdd_
//void zgesdd_

// Eigen values/vectors
//void dggev_
//void zggev_
]]

-- complex versions

return lapack