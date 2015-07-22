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

ffi.cdef[[
// Y = alpha * A * X + beta * Y
void dgemv_(char* TRANS, const int* M, const int* N,
           double* alpha, double* A, const int* LDA, double* X,
           const int* INCX, double* beta, double* C, const int* INCY);
// C = alpha * A * B + beta * C                                               
void dgemm_(char* TRANSA, char* TRANSB, const int* M,
           const int* N, const int* K, double* alpha, double* A,
           const int* LDA, double* B, const int* LDB, double* beta,
           double* C, const int* LDC);

// N x N
void dgesv_(const int *N, const int *nrhs, double *A, const int *lda,
            int *ipiv, double *b, const int *ldb, int *info);
//void zgesv_

// Driver N x N
//void dgesvx_
//void zgesvx_

// M x N
void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
    		double *A, const int *lda, double *b, const int *ldb, double *work,
    		const int *lwork, int *info);
//void zgels_

// Driver M x N: min | b - Ax | using SVD
//void dgelss_
//void zgelss_


// SVD
void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N,
        	 double* A, const int* lda, double* S, double* U, const int* ldu,
        	 double* VT, const int* ldvt, double* work,const int* lwork,
        	 const int* info);
//void zgesvd_

// Eigen values/vectors
//void dggev_
//void zggev_
]]

return lapack