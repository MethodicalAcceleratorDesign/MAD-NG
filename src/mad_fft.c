/*
 o-----------------------------------------------------------------------------o
 |
 | FFT module implementation
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
//#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"

// --- implementation ---------------------------------------------------------o

#define CHKX   assert( x )
#define CHKXR  assert( x && r )

// -- FFT ---------------------------------------------------------------------o

#include <fftw3.h>

void // x [n] -> r [n]
mad_vec_fft (const num_t x[], cnum_t r[], ssz_t n)
{
  CHKXR;
  mad_alloc_tmp(cnum_t, cx, n);
  mad_vec_copyv(x, cx, n);
  mad_cvec_fft(cx, r, n);
  mad_free_tmp(cx);
}

void // x [n] -> r [n/2+1]
mad_vec_rfft (const num_t x[], cnum_t r[], ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_r2c_1d(n, (num_t*)x, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cvec_fft (const cnum_t x[], cnum_t r[], ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_1d(n, (cnum_t*)x, r, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cvec_ifft(const cnum_t x[], cnum_t r[], ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_1d(n, (cnum_t*)x, r, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_cvec_muln(r, 1.0/n, r, n);
}

void // x [n/2+1] -> r [n]
mad_cvec_irfft (const cnum_t x[], num_t r[], ssz_t n)
{
  CHKXR;
  ssz_t nn = n/2+1;
  mad_alloc_tmp(cnum_t, cx, nn);
  mad_cvec_copy(x, cx, nn);
  fftw_plan p = fftw_plan_dft_c2r_1d(n, cx, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_free_tmp(cx);
  mad_vec_muln(r, 1.0/n, r, n);
}

void // x [m x n] -> r [m, n/2+1]
mad_mat_fft (const num_t x[], cnum_t r[], ssz_t m, ssz_t n)
{
  CHKXR;
  mad_alloc_tmp(cnum_t, cx, m*n);
  mad_vec_copyv(x, cx, m*n);
  mad_cmat_fft(cx, r, m, n);
  mad_free_tmp(cx);
}

void // x [m x n] -> r [m, n/2+1]
mad_mat_rfft (const num_t x[], cnum_t r[], ssz_t m, ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_r2c_2d(m, n, (num_t*)x, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cmat_fft (const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_2d(m, n, (cnum_t*)x, r, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cmat_ifft(const cnum_t x[], cnum_t r[], ssz_t m, ssz_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_2d(m, n, (cnum_t*)x, r, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_cvec_muln(r, 1.0/(m*n), r, m*n);
}

void // x [m x n/2+1] -> r [m x n]
mad_cmat_irfft (const cnum_t x[], num_t r[], ssz_t m, ssz_t n)
{
  CHKXR;
  ssz_t nn = m*(n/2+1);
  mad_alloc_tmp(cnum_t, cx, nn);
  mad_cvec_copy(x, cx, nn);
  fftw_plan p = fftw_plan_dft_c2r_2d(m, n, cx, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_free_tmp(cx);
  mad_vec_muln(r, 1.0/(m*n), r, m*n);
}

/* -- NFFT --------------------------------------------------------------------o
[1] J. Keiner, S. Kunis and D. Potts, "Using NFFT 3 — A Software Library for
    Various Nonequispaced Fast Fourier Transforms", ACM Transactions on
    Mathematical Software, Vol. 36, No. 4, Article 19, Pub. date: August 2009.
[2] J. Keiner, S. Kunis and D. Potts, "NFFT 3.0 - Tutorial",
    http://www.tu-chemnitz.de/~potts/nfft
[3] S. Kunis and D. Potts, "Time and memory requirements of the Nonequispaced
    FFT", report, 2006.
[4] L. Kammerer, S. Kunis, I. Melzer, D. Potts, and T. Volkmer, "Computational
    Methods for the Fourier Analysis of Sparse High-Dimensional Functions",
    Springer Lecture Notes 102, 2014.
[5] http://github.com/NFFT/nfft

Notations vs [1]: X=f_hat, Y=f
eq. 2.1 (forward DFT, i.e. time to frequency):
Y_j = sum_{k=0}^{N-1} X_k e^{-2\pi i k j/N},  (j=0,…,N-1)

eq. 2.2 (backward DFT, frequency to time):
X_k = sum_{j=0}^{N-1} Y_j e^{+2\pi i k j/N},  (k=0,…,N-1)

eq. 2.3 (forward NDFT, i.e. time to frequency):
Y_j = sum_{k=0}^{N-1} X_k e^{-2\pi i k x_j},  (j=0,…,M-1)
if M = N and x_j = j/N, eq. 2.3 becomes eq. 2.1

eq. 2.4 (adjoint NDFT, i.e. frequency to time):
X_k = sum_{j=0}^{M-1} Y_j e^{+2\pi i k j/N},  (k=0,…,N-1)
if M = N and x_j = j/N, eq. 2.4 becomes eq. 2.2

Expected use:
eq. 2.3 (forward NDFT, i.e. time to frequency):
Y_k = sum_{j=0}^{M-1} X_j e^{-2\pi i k x_j},  (k=0,…,N-1)
if M = N and x_j = j/N, eq. 2.3 becomes eq. 2.1

eq. 2.4 (adjoint NDFT, i.e. frequency to time):
X_j = sum_{k=0}^{N-1} Y_k e^{+2\pi i k x_j},  (j=0,…,M-1)
if M = N and x_j = j/N, eq. 2.4 becomes eq. 2.2

Application of NFFT3 (with FFTW3 compatibility):
forward NFFT, i.e. time to frequency:
  - Use nfft_adjoint with negated input time nodes and shifted output frequency
    indexes by length/2.
backward NFFT, i.e. frequency to time:
  - Use nfft_trafo with negated input time nodes and shifted input frequency
    indexes by length/2, and normalized output signal by 1/length.
*/

#ifndef MAD_NO_NFFT
#include <nfft3.h>

static nfft_plan p;
static ssz_t p_n, p_n1, p_n2, p_m;

void
mad_vec_nfft (const num_t x[], const num_t x_node[], cnum_t r[], ssz_t n, ssz_t nr)
{
  CHKX;
  mad_alloc_tmp(cnum_t, cx, n);
  mad_vec_copyv(x, cx, n);
  mad_cvec_nfft(cx, x_node, r, n, nr);
  mad_free_tmp(cx);
}

void // time to frequency
mad_cvec_nfft (const cnum_t x[], const num_t x_node[], cnum_t r[], ssz_t n, ssz_t nr)
{
  assert( x && r );
  int precomp = 0;
  if (n != p_n || nr != p_m) {
    nfft_finalize(&p);
    nfft_init_1d (&p, n, nr);
    p_n = n, p_m = nr, precomp = 1;
  }
  if (x_node || precomp) {
    for (idx_t i=0; i < n; i++)  // adjoint transform needs -x_node
      p.x[i] = x_node[i] == -0.5 ? 0.4999999999999999 : -x_node[i];
    if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  }
  mad_cvec_copy(x, p.f, n);
  const char *error_str = nfft_check(&p);
  if (error_str) error("%s", error_str);
  nfft_adjoint(&p); // nfft_adjoint_direct(&p);
  mad_cvec_copy(p.f_hat+nr/2, r, nr/2); // for compatibility with FFTW
  mad_cvec_copy(p.f_hat, r+nr/2, nr/2);
}

void // frequency to time
mad_cvec_infft (const cnum_t x[], const num_t r_node[], cnum_t r[], ssz_t n, ssz_t nx)
{
  assert( x && r );
  int precomp = 0;
  if (n != p_n || nx != p_m) {
    nfft_finalize(&p);
    nfft_init_1d (&p, n, nx);
    p_n = n, p_m = nx, precomp = 1;
  }
  if (r_node || precomp) {
    for (idx_t i=0; i < n; i++) // forward transform needs -r_node
      p.x[i] = r_node[i] == -0.5 ? 0.4999999999999999 : -r_node[i];
    if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  }
  mad_cvec_copy(x+nx/2, p.f_hat, nx/2); // for compatibility with FFTW
  mad_cvec_copy(x, p.f_hat+nx/2, nx/2);
  const char *error_str = nfft_check(&p);
  if (error_str) error("%s", error_str);
  nfft_trafo(&p); // nfft_trafo_direct(&p);
  mad_cvec_copy(p.f, r, n);
  mad_cvec_muln(r, 1.0/n, r, n);
}

void
mad_mat_nfft (const num_t x[], const num_t x_node[], cnum_t r[], ssz_t m, ssz_t n, ssz_t nr)
{
  CHKX;
  mad_alloc_tmp(cnum_t, cx, m*n);
  mad_vec_copyv(x, cx, m*n);
  mad_cmat_nfft(cx, x_node, r, m, n, nr);
  mad_free_tmp(cx);
}

void // space to frequency
mad_cmat_nfft (const cnum_t x[], const num_t x_node[], cnum_t r[], ssz_t m, ssz_t n, ssz_t nr)
{
  assert( x && r );
  int precomp = 0;
  if (m != p_n1 || n != p_n2 || nr != p_m) {
    nfft_finalize(&p);
    nfft_init_2d (&p, m, n, nr);
    p_n1 = m, p_n2 = n, p_m = nr, precomp = 1;
  }
  if (x_node || precomp) {
    for (ssz_t i=0; i < m*n; i++)  // adjoint transform needs -x_node
      p.x[i] = x_node[i] == -0.5 ? 0.4999999999999999 : -x_node[i];
    if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  }
  mad_cvec_copy(x, p.f, m*n);
  const char *error_str = nfft_check(&p);
  if (error_str) error("%s", error_str);
  nfft_adjoint(&p); // nfft_adjoint_direct(&p);
  mad_cvec_copy(p.f_hat+nr/2, r, nr/2); // for compatibility with FFTW ?? (TBC)
  mad_cvec_copy(p.f_hat, r+nr/2, nr/2);
}

void // frequency to space
mad_cmat_infft (const cnum_t x[], const num_t r_node[], cnum_t r[], ssz_t m, ssz_t n, ssz_t nx)
{
  assert( x && r );
  int precomp = 0;
  if (m != p_n1 || n != p_n2 || nx != p_m) {
    nfft_finalize(&p);
    nfft_init_2d (&p, m, n, nx);
    p_n1 = m, p_n2 = n, p_m = nx, precomp = 1;
  }
  if (r_node || precomp) {
    for (ssz_t i=0; i < m*n; i++) // forward transform needs -r_node
      p.x[i] = r_node[i] == -0.5 ? 0.4999999999999999 : -r_node[i];
    if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  }
  mad_cvec_copy(x+nx/2, p.f_hat, nx/2); // for compatibility with FFTW ?? (TBC)
  mad_cvec_copy(x, p.f_hat+nx/2, nx/2);
  const char *error_str = nfft_check(&p);
  if (error_str) error("%s", error_str);
  nfft_trafo(&p); // nfft_trafo_direct(&p);
  mad_cvec_copy(p.f, r, m*n);
  mad_cvec_muln(r, 1.0/(m*n), r, m*n);
}

#endif // MAD_NO_NFFT

// -- CLEANUP -----------------------------------------------------------------o

void
mad_fft_cleanup(void)
{
#ifndef MAD_NO_NFFT
  nfft_finalize(&p);  memset(&p, 0, sizeof p);
#endif
  fftw_cleanup();
}
