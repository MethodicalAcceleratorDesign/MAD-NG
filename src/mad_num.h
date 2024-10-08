#ifndef MAD_NUM_H
#define MAD_NUM_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Number module interface
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

  Purpose:
  - wrappers around functions of real and complex numbers for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include <math.h>
#include "mad_def.h"

// --- interface --------------------------------------------------------------o

int   mad_num_sign     (num_t x); // -1, 0, 1
int   mad_num_sign1    (num_t x); // -1, 1

num_t mad_num_fact     (int n);   // n in Z -> n!
num_t mad_num_fact2    (int n);   // n in Z -> n!! (wiki/Double_factorial)
num_t mad_num_binom    (int n, int k); // (n,k) in Z^2 -> n!/(k!(n-k)!)

static inline
num_t mad_num_div      (num_t x, num_t y) { return x/y; }
static inline
num_t mad_num_inv      (num_t x)          { return 1/x; }
num_t mad_num_sinc     (num_t x);
num_t mad_num_sinhc    (num_t x);
num_t mad_num_asinc    (num_t x);
num_t mad_num_asinhc   (num_t x);
num_t mad_num_powi     (num_t x, int n);

cpx_t mad_cpx_div      (cpx_t x, cpx_t y);
cpx_t mad_cpx_inv      (cpx_t x);
cpx_t mad_cpx_sinc     (cpx_t x);
cpx_t mad_cpx_sinhc    (cpx_t x);
cpx_t mad_cpx_asinc    (cpx_t x);
cpx_t mad_cpx_asinhc   (cpx_t x);
cpx_t mad_cpx_powi     (cpx_t x, int n);

num_t mad_cpx_abs_r    (num_t x_re, num_t x_im);
num_t mad_cpx_arg_r    (num_t x_re, num_t x_im);

void  mad_cpx_sqrt_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_exp_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_log_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_log10_r  (num_t x_re, num_t x_im, cpx_t *r);

void  mad_cpx_sin_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_cos_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_tan_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_sinh_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_cosh_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_tanh_r   (num_t x_re, num_t x_im, cpx_t *r);

void  mad_cpx_asin_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_acos_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_atan_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_asinh_r  (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_acosh_r  (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_atanh_r  (num_t x_re, num_t x_im, cpx_t *r);

void  mad_cpx_sinc_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_sinhc_r  (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_asinc_r  (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_asinhc_r (num_t x_re, num_t x_im, cpx_t *r);

void  mad_cpx_unit_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_proj_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_rect_r   (num_t  rho, num_t  ang, cpx_t *r);
void  mad_cpx_polar_r  (num_t x_re, num_t x_im, cpx_t *r);

void  mad_cpx_invsqrt_r(num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_inv_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_div_r    (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r);
void  mad_cpx_mod_r    (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r);
void  mad_cpx_pow_r    (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r);
void  mad_cpx_powi_r   (num_t x_re, num_t x_im, int   n,                cpx_t *r);

// --- Faddeeva based functions -----------------------------------------------o

num_t mad_num_wf       (num_t x);
num_t mad_num_erf      (num_t x);
num_t mad_num_erfc     (num_t x);
num_t mad_num_erfi     (num_t x);
num_t mad_num_erfcx    (num_t x);
num_t mad_num_dawson   (num_t x);

cpx_t mad_cpx_wf       (cpx_t x);
cpx_t mad_cpx_erf      (cpx_t x);
cpx_t mad_cpx_erfc     (cpx_t x);
cpx_t mad_cpx_erfi     (cpx_t x);
cpx_t mad_cpx_erfcx    (cpx_t x);
cpx_t mad_cpx_dawson   (cpx_t x);

void  mad_cpx_wf_r     (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_erf_r    (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_erfc_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_erfi_r   (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_erfcx_r  (num_t x_re, num_t x_im, cpx_t *r);
void  mad_cpx_dawson_r (num_t x_re, num_t x_im, cpx_t *r);

// --- RNG --------------------------------------------------------------------o

// MAD peudo-random number generator

typedef struct prng_state_ prng_state_t;

num_t mad_num_rand     (prng_state_t*);             // [0.,1.)
u64_t mad_num_randi    (prng_state_t*);             // [0,ULLONG_MAX]
void  mad_num_randjump (prng_state_t*);
void  mad_num_randseed (prng_state_t*, num_t seed);

// --- private
struct prng_state_ {
  u64_t s[4];
};

// MAD-X pseudo-random number generator
typedef struct xrng_state_ xrng_state_t;

num_t mad_num_xrand     (xrng_state_t*);           // [0.,1.)
u32_t mad_num_xrandi    (xrng_state_t*);           // [0,UINT_MAX]
void  mad_num_xrandseed (xrng_state_t*, u32_t seed);

// --- private
struct xrng_state_ {
  int s[55];
  idx_t n;
};

// --- OMP --------------------------------------------------------------------o

num_t mad_num_suminv (u64_t n); // dummy function for testing OpenMP

// ----------------------------------------------------------------------------o

#endif // MAD_NUM_H
