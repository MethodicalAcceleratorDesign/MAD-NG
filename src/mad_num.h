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
#include "mad_defs.h"

// --- interface --------------------------------------------------------------o

int    mad_num_sign    (num_t x); // -1, 0, 1
int    mad_num_sign1   (num_t x); // -1, 1

num_t  mad_num_sinc    (num_t x);
num_t  mad_num_sinhc   (num_t x);
num_t  mad_num_powi    (num_t x, int n);

cnum_t mad_cnum_sinc   (cnum_t x);
cnum_t mad_cnum_sinhc  (cnum_t x);
cnum_t mad_cnum_powi   (cnum_t x, int n);

num_t mad_cnum_abs_r   (num_t x_re, num_t x_im);
num_t mad_cnum_arg_r   (num_t x_re, num_t x_im);

void  mad_cnum_sqrt_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_exp_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_log_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_log10_r (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_sin_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_cos_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_tan_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_sinh_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_cosh_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_tanh_r  (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_asin_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_acos_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_atan_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_asinh_r (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_acosh_r (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_atanh_r (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_sinc_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_sinhc_r (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_unit_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_proj_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_rect_r  (num_t  rho, num_t  ang, cnum_t *r);
void  mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_inv_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_div_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_mod_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_pow_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_powi_r  (num_t x_re, num_t x_im, int   n,                cnum_t *r);

// --- Faddeeva based functions -----------------------------------------------o

num_t   mad_num_erf    (num_t x);
num_t   mad_num_erfc   (num_t x);
num_t   mad_num_erfi   (num_t x);
num_t   mad_num_erfw   (num_t x);
num_t   mad_num_erfcx  (num_t x);
num_t   mad_num_dawson (num_t x);

cnum_t  mad_cnum_erf   (cnum_t x, num_t relerr);
cnum_t  mad_cnum_erfc  (cnum_t x, num_t relerr);
cnum_t  mad_cnum_erfi  (cnum_t x, num_t relerr);
cnum_t  mad_cnum_erfw  (cnum_t x, num_t relerr);
cnum_t  mad_cnum_erfcx (cnum_t x, num_t relerr);
cnum_t  mad_cnum_dawson(cnum_t x, num_t relerr);

void  mad_cnum_erf_r   (num_t x_re, num_t x_im, num_t relerr, cnum_t *r);
void  mad_cnum_erfc_r  (num_t x_re, num_t x_im, num_t relerr, cnum_t *r);
void  mad_cnum_erfi_r  (num_t x_re, num_t x_im, num_t relerr, cnum_t *r);
void  mad_cnum_erfw_r  (num_t x_re, num_t x_im, num_t relerr, cnum_t *r);
void  mad_cnum_erfcx_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r);
void  mad_cnum_dawson_r(num_t x_re, num_t x_im, num_t relerr, cnum_t *r);

// --- RNG --------------------------------------------------------------------o

typedef struct rng_state rng_state_t;

num_t mad_num_rand     (rng_state_t*);             // [0.,1.)
u64_t mad_num_randi    (rng_state_t*);             // [0,ULLONG_MAX]
void  mad_num_randseed (rng_state_t*, num_t seed);
void  mad_num_randjump (rng_state_t*, const rng_state_t* ref);

// --- private
struct rng_state {
  u64_t s[16];
  int p;
};

// ----------------------------------------------------------------------------o

#endif // MAD_NUM_H
