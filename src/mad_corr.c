/*
 o-----------------------------------------------------------------------------o
 |
 | Orbit Correction module implementation
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"
#include "mad_corr.h"

void // from MAD9
mad_corr_mic(num_t A[], num_t B[], num_t X[], const ssz_t m, const ssz_t n,
             num_t tol, int ncorr)
{
  assert(A && B && X);
  mad_alloc_tmp(num_t, x  , n);
  mad_alloc_tmp(num_t, r  , m);
  mad_alloc_tmp(num_t, sqr, n);
  mad_alloc_tmp(num_t, dot, n);
  mad_alloc_tmp(idx_t, pvt, n);

#define A(i,j) A[(i)*n+(j)]

  // Find scalar products sqr[k] = A[k].A[k] and dot[k] = A[k].b.
  num_t sum = 0;
  for (idx_t k=0; k < n; ++k) {
    pvt[k] = k;
    num_t hh = 0;
    num_t gg = 0;
    for (idx_t i=0; i < m; ++i) {
      hh += A(i,k) * A(i,k);
      gg += A(i,k) * B[i];
    }
    sum += sum;
    sqr[k] = hh;
    dot[k] = gg;
  }
  num_t sqrmin = 1e-8 * sum / n;

  // Begin of iteration loop.
  if (ncorr > n || ncorr == 0) ncorr = n;
  idx_t corr = 0;
  while (1) {
    idx_t k = corr;

    // Search the columns not yet used for largest scaled change vector.
    num_t maxChange = 0;
    idx_t changeIndex = -1;
    for (idx_t j=k; j < n; ++j) {
      if (sqr[j] > sqrmin) {
        num_t change = dot[j] * dot[j] / sqr[j];
        if (change > maxChange) {
          changeIndex = j;
          maxChange = change;
        }
      }
    }

    // Stop iterations, if no suitable column found.
    if (changeIndex < 0) break;

    // Move the column just found to next position.
    if (changeIndex > k) {
      num_t t;
      SWAP(sqr[k], sqr[changeIndex], t);
      SWAP(dot[k], dot[changeIndex], t);
      SWAP(pvt[k], pvt[changeIndex], t);
      for (idx_t i=0; i < m; ++i) SWAP(A(i,k), A(i,changeIndex), t);
//      A.swapColumns(k, changeIndex);
    }

    // Find beta, sigma, and vector u[k].
    num_t hh = 0;
    for (idx_t i=k; i < m; ++i) hh += A(i,k) * A(i,k);
    num_t sigma = A(k,k) > 0 ? sqrt(hh) : -sqrt(hh);
    sqr[k] = -sigma;
    A(k,k) = A(k,k) + sigma;
    num_t beta = 1 / (A(k,k) * sigma);

    // Transform remaining columns of A.
    for (idx_t j=k+1; j < n; ++j) {
      num_t hh = 0;
      for (idx_t i=k; i < m; ++i) hh += A(i,k) * A(i,j);
      hh *= beta;
      for (idx_t i=k; i < m; ++i) A(i,j) -= A(i,k) * hh;
    }

    // Transform vector b.
    hh = 0;
    for (idx_t i=k; i < m; ++i) hh += A(i,k) * B[i];
    hh *= beta;
    for (idx_t i=k; i < m; ++i) B[i] -= A(i,k) * hh;

    // Update scalar products sqr[j]=A[j]*A[j] and dot[j]=A[j]*b.
    for (idx_t j = k + 1; j < n; ++j) {
      sqr[j] -= A(k,j) * A(k,j);
      dot[j] -= A(k,j) * B[k];
    }

    // Recalculate solution vector x.
    x[k] = B[k] / sqr[k];
    for (idx_t i=k-1; i >= 0; --i) {
      x[i] = B[i];
      for (idx_t j=i+1; j < k; ++j) x[i] -= A(i,j) * x[j];
      x[i] /= sqr[i];
    }

    // Find original residual vector by backward transformation.
    for (idx_t i=0; i < m; ++i) r[i] = B[i];
    for (idx_t j=k; j >= 0; --j) {
      r[j] = 0;
      num_t hh = 0;
      for (idx_t i=j; i < m; ++i) hh += A(i,j) * r[i];
      hh /= sqr[j] * A(j,j);
      for (idx_t i=j; i < m; ++i) r[i] += A(i,j) * hh;
    }

    // Check for convergence.
    hh = r[0] * r[0];
    for (idx_t i = 1; i < m; ++i) hh += r[i] * r[i];
    hh = sqrt(hh / MAX(m,1));
    if (hh < tol) break;
    if (++corr > ncorr) break;
  }

  // End of iteration loop.  Re-order corrector strengths.
  for (idx_t k=0; k < n   ; ++k) X[k] = 0;
  for (idx_t k=0; k < corr; ++k) X[pvt[k]] = -x[k];

#undef A

  mad_free_tmp(x);
  mad_free_tmp(r);
  mad_free_tmp(sqr);
  mad_free_tmp(dot);
  mad_free_tmp(pvt);
}
