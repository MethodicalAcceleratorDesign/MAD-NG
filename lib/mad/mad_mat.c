// gcc -std=c11 -W -Wall -Wextra -pedantic -O3 -ffast-math -ftree-vectorize -shared -fPIC -static-libgcc *.c -o libmad-OSX.so -lm

#include "mad_mat.h"
#include <complex.h>
#include <assert.h>

typedef double           num_t;
typedef double _Complex cnum_t;

#define CHKXR  assert( x && y && r && x != r )
#define CHKYR  assert( x && y && r && y != r )
#define CHKXYR assert( x && y && r && x != r && y != r )

// [m x n] = [m x p] * [p x n]
/*
#define MMUL \
  for (size_t i=0; i < m; i++) \
    for (size_t j=0; j < n; j++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
*/

#define MMUL \
  for (size_t i=0; i < m; i++) { \
    for (size_t j=0; j < n-3; j+=4) { \
      r[i*n+j  ] = 0; \
      r[i*n+j+1] = 0; \
      r[i*n+j+2] = 0; \
      r[i*n+j+3] = 0; \
      for (size_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
        r[i*n+j+2] += x[i*p+k] * y[k*n+j+2]; \
        r[i*n+j+3] += x[i*p+k] * y[k*n+j+3]; \
      } \
    } \
  } \
  if (n & 2) { \
    size_t j=n-3; \
    for (size_t i=0; i < m; i++) { \
      r[i*n+j  ] = 0; \
      r[i*n+j+1] = 0; \
      for (size_t k=0; k < p; k++) { \
        r[i*n+j  ] += x[i*p+k] * y[k*n+j  ]; \
        r[i*n+j+1] += x[i*p+k] * y[k*n+j+1]; \
      } \
    } \
  } \
  if (n & 1) { \
    size_t j=n-1; \
    for (size_t i=0; i < m; i++) { \
      r[i*n+j] = 0; \
      for (size_t k=0; k < p; k++) \
        r[i*n+j] += x[i*p+k] * y[k*n+j]; \
    } \
  } \

// [m x 1] = [m x p] * [p x 1]
#define MULV \
  for (size_t i=0; i < m; i++) { \
    r[i] = 0; \
    for (size_t k=0; k < p; k++) \
      r[i] += x[i*p+k] * y[k]; \
  } \

// [1 x n] = [1 x p] * [p x n]
#define VMUL \
  for (size_t j=0; j < n; j++) { \
    r[j] = 0; \
    for (size_t k=0; k < p; k++) \
      r[j] += x[k] * y[k*n+j]; \
  } \


void mad_mat_mul (const num_t *x, const num_t *y, num_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; MMUL; }

void mad_mat_mulm (const num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKYR; MMUL; }

void mad_mat_mulv (const num_t *x, const num_t *y, num_t *r, size_t m, size_t p)
{ CHKYR; MULV; }

void mad_mat_mulc (const num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p)
{ CHKYR; MULV; }

void mad_mat_vmul (const num_t *x, const num_t *y, num_t *r, size_t n, size_t p)
{ CHKXR; VMUL; }

void mad_mat_cmul (const cnum_t *x, const num_t *y, cnum_t *r, size_t n, size_t p)
{ CHKXR; VMUL; }

void mad_cmat_mul (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXYR; MMUL; }

void mad_cmat_mulm (const cnum_t *x, const num_t *y, cnum_t *r, size_t m, size_t n, size_t p)
{ CHKXR; MMUL; }

void mad_cmat_mulv (const cnum_t *x, const num_t *y, cnum_t *r, size_t m, size_t p)
{ CHKXR; MULV; }

void mad_cmat_mulc (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p)
{ CHKYR; MULV; }

void mad_cmat_vmul (const num_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p)
{ CHKYR; VMUL; }

void mad_cmat_cmul (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p)
{ CHKXR; VMUL; }

void mad_mat_trans (const num_t *x, num_t *r, size_t m, size_t n)
{ 
  assert(x && r);
  if (x == r) {
    assert(m == n);
    for (size_t i=0  ; i < m; i++)
    for (size_t j=i+1; j < n; j++)
      r[j*n+i] = r[i*n+j];
  } else {
    for (size_t i=0; i < m; i++)
    for (size_t j=0; j < n; j++)
      r[j*m+i] = x[i*n+j];
  }
}

void mad_cmat_trans (const cnum_t *x, cnum_t *r, size_t m, size_t n)
{ 
  assert(x && r);
  if (x == r) {
    assert(m == n);
    for (size_t i=0; i < m; i++)
    for (size_t j=i; j < n; j++)
      r[j*n+i] = conj(r[i*n+j]);
  } else {
    for (size_t i=0; i < m; i++)
    for (size_t j=0; j < n; j++)
      r[j*m+i] = conj(x[i*n+j]);
  }
}
