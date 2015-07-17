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
// C = alpha * A * B + beta * C                                               
void dgemm_(char* TRANSA, char* TRANSB, const int* M,
           const int* N, const int* K, double* alpha, double* A,
           const int* LDA, double* B, const int* LDB, double* beta,
           double* C, const int* LDC);
// Y = alpha * A * X + beta * Y                                               
void dgemv_(char* TRANS, const int* M, const int* N,
           double* alpha, double* A, const int* LDA, double* X,
           const int* INCX, double* beta, double* C, const int* INCY);
]]

return lapack