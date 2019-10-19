/*
 o-----------------------------------------------------------------------------o
 |
 | Number module implementation
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
#include <complex.h>
#include <assert.h>

#include "mad_num.h"
#include "mad_cst.h"

// --- implementation ---------------------------------------------------------o

#define CHKR  assert( r )

#define CNUM2(a,b) (* (cnum_t*) & (num_t[2]) { a, b })
#define CNUM(a) CNUM2(MKNAME(a,_re), MKNAME(a,_im))

// --- num

int mad_num_sign (num_t x)
{
  return (x > 0) - (x < 0);  // -1, 0, 1
}

int mad_num_sign1 (num_t x)
{
  return ((int[]){ -1, 1 })[!signbit(x)]; // -1, 1: works for ±0, ±inf and ±NaN
}

num_t mad_num_sinc (num_t x)
{
  return fabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x : sin (x)/x;
}

num_t mad_num_sinhc (num_t x)
{
  return fabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x : sinh(x)/x;
}

num_t mad_num_powi (num_t x, int n)
{
  num_t r = 1;
  if (n < 0) n = -n, x = 1/x;
  for (;;) {
    if (n &   1) r *= x;
    if (n >>= 1) x *= x; else break;
  }
  return r;
}

// --- cnum

cnum_t mad_cnum_sinc  (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x : csin (x)/x;
}

cnum_t mad_cnum_sinhc (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x : csinh(x)/x;
}

cnum_t mad_cnum_powi (cnum_t x, int n)
{
  cnum_t r = 1;
  if (n < 0) n = -n, x = 1/x;
  for (;;) {
    if (n &   1) r *= x;
    if (n >>= 1) x *= x; else break;
  }
  return r;
}

num_t mad_cnum_abs_r  (num_t x_re, num_t x_im) { return cabs( CNUM(x) ); }
num_t mad_cnum_arg_r  (num_t x_re, num_t x_im) { return carg( CNUM(x) ); }

void mad_cnum_sqrt_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = csqrt  ( CNUM(x) ); }
void mad_cnum_exp_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cexp   ( CNUM(x) ); }
void mad_cnum_log_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = clog   ( CNUM(x) ); }
void mad_cnum_log10_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = clog   ( CNUM(x) )/log(10); }

void mad_cnum_sin_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = csin   ( CNUM(x) ); }
void mad_cnum_cos_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = ccos   ( CNUM(x) ); }
void mad_cnum_tan_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = ctan   ( CNUM(x) ); }
void mad_cnum_sinh_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = csinh  ( CNUM(x) ); }
void mad_cnum_cosh_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = ccosh  ( CNUM(x) ); }
void mad_cnum_tanh_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = ctanh  ( CNUM(x) ); }

void mad_cnum_asin_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = casin  ( CNUM(x) ); }
void mad_cnum_acos_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cacos  ( CNUM(x) ); }
void mad_cnum_atan_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = catan  ( CNUM(x) ); }
void mad_cnum_asinh_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = casinh ( CNUM(x) ); }
void mad_cnum_acosh_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cacosh ( CNUM(x) ); }
void mad_cnum_atanh_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = catanh ( CNUM(x) ); }

void mad_cnum_proj_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cproj  ( CNUM(x) ); }

void mad_cnum_sinc_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_sinc ( CNUM(x) ); }
void mad_cnum_sinhc_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_sinhc( CNUM(x) ); }

void mad_cnum_powi_r  (num_t x_re, num_t x_im, int n, cnum_t *r)
{ CHKR; *r = mad_cnum_powi( CNUM(x), n ); }

void mad_cnum_unit_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / cabs( CNUM(x) ); }

void mad_cnum_rect_r (num_t rho, num_t ang, cnum_t *r)
{ CHKR; *r = CNUM2( rho * cos(ang), rho * sin(ang) ); }

void mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM2( cabs(CNUM(x)), carg(CNUM(x)) ); }

void mad_cnum_inv_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = 1 / CNUM(x); }

void mad_cnum_div_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / CNUM(y);  }

void mad_cnum_mod_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; cnum_t cr = CNUM(x) / CNUM(y);
  *r = CNUM(x) - CNUM(y) * CNUM2(round(creal(cr)), round(cimag(cr))); }

void mad_cnum_pow_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = cpow( CNUM(x), CNUM(y) ); }

// --- Faddeeva function and variants from MIT --------------------------------o

#include "mad_erfw.h"

num_t  mad_num_erf    (num_t x) { return Faddeeva_erf_re   (x); }
num_t  mad_num_erfc   (num_t x) { return Faddeeva_erfc_re  (x); }
num_t  mad_num_erfi   (num_t x) { return Faddeeva_erfi_re  (x); }
num_t  mad_num_erfcx  (num_t x) { return Faddeeva_erfcx_re (x); }
num_t  mad_num_dawson (num_t x) { return Faddeeva_Dawson_re(x); }
num_t  mad_num_erfw   (num_t x) { return Faddeeva_w_im     (x); }

cnum_t mad_cnum_erf   (cnum_t x, num_t relerr) { return Faddeeva_erf   (x, relerr); }
cnum_t mad_cnum_erfc  (cnum_t x, num_t relerr) { return Faddeeva_erfc  (x, relerr); }
cnum_t mad_cnum_erfi  (cnum_t x, num_t relerr) { return Faddeeva_erfi  (x, relerr); }
cnum_t mad_cnum_erfcx (cnum_t x, num_t relerr) { return Faddeeva_erfcx (x, relerr); }
cnum_t mad_cnum_dawson(cnum_t x, num_t relerr) { return Faddeeva_Dawson(x, relerr); }
cnum_t mad_cnum_erfw  (cnum_t x, num_t relerr) { return Faddeeva_w     (x, relerr); }

void mad_cnum_erf_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erf (CNUM(x), relerr); }

void mad_cnum_erfc_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfc (CNUM(x), relerr); }

void mad_cnum_erfi_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfi (CNUM(x), relerr); }

void mad_cnum_erfcx_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfcx (CNUM(x), relerr); }

void mad_cnum_dawson_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_Dawson (CNUM(x), relerr); }

void mad_cnum_erfw_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_w (CNUM(x), relerr); }

// -- RNG XorShift1024* -------------------------------------------------------o

enum {
  N = 16,
  rand_state_expected_size = N * sizeof(u64_t),
  rand_state_actual_size   = offsetof(rng_state_t,p),
  invalid_rand_state_size  = 1/(rand_state_actual_size == rand_state_expected_size)
};

union numbit {
  u64_t u;
  num_t d;
};

u64_t mad_num_randi (rng_state_t *restrict st)
{
  int p = st->p;
  u64_t *s = st->s;
  const u64_t s0 = s[p];
  u64_t s1 = s[p = (p + 1) & (N - 1)];
  s1 ^= s1 << 31; // A
  s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // B, C
  st->p = p;
  return s[p] * 1181783497276652981ULL; // number within [0,ULLONG_MAX]
}

num_t mad_num_rand (rng_state_t *restrict st)
{
  u64_t r = mad_num_randi(st);
  r = (r & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL;
  const union numbit n = { .u = r }; // number within [1.,2.)
  return n.d - 1.0;                  // number within [0.,1.)
}

void mad_num_randseed (rng_state_t *restrict st, num_t seed)
{
  u64_t *s = st->s;
  const union numbit n = { .d = seed ? seed : M_PI }; // avoid zero
  s[0] = n.u * 33;
  for (int i = 1; i < N; i++)
    s[i] = s[i-1] * 33;
  for (int i = 0; i < N; i++)
    mad_num_randi(st);

  mad_num_randjump(st,0);
}

void mad_num_randjump (      rng_state_t *restrict st,
                       const rng_state_t *restrict ref)
{
  static const u64_t jump[N] = {
    0x84242f96eca9c41dULL, 0xa3c65b8776f96855ULL, 0x5b34a39f070b5837ULL,
    0x4489affce4f31a1eULL, 0x2ffeeb0a48316f40ULL, 0xdc2d9891fe68c022ULL,
    0x3659132bb12fea70ULL, 0xaac17d8efa43cab8ULL, 0xc4cb815590989b13ULL,
    0x5ee975283d71c93bULL, 0x691548c86c1bd540ULL, 0x7910c41d10a1e6a5ULL,
    0x0b5fc64563b3e2a8ULL, 0x047f7684e9fc949dULL, 0xb99181f2d8f685caULL,
    0x284600e3f30e38c3ULL
  };

  if (ref)
    for (int i = 0; i < N; i++) {
      st->s[i] = ref->s[i];
      st->p    = ref->p;
    }

  int p = st->p;
  u64_t *s = st->s;
  u64_t t[N] = { 0 };

  for(int i = 0; i < N; i++)
    for(int b = 0; b < 64; b++) {
      if (jump[i] & 1ULL << b)
        for(int j = 0; j < N; j++)
          t[j] ^= s[(j + p) & (N - 1)];
      mad_num_randi(st);
    }
  for(int j = 0; j < N; j++)
    s[(j + p) & (N - 1)] = t[j];
}

#undef N
