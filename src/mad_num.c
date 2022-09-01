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
#include <float.h>
#include <complex.h>
#include <stdlib.h>
#include <assert.h>

#include "mad_num.h"
#include "mad_cst.h"

// --- helpers ----------------------------------------------------------------o

static inline num_t
fact(int n)
{
  static num_t f[171] = {1, 0};

  if (!f[1])
    for (int i=1; i < 171; ++i) f[i] = i*f[i-1];

  return n < 171 ? f[n] : INFINITY;
}

// --- implementation ---------------------------------------------------------o

#define CHKR  assert( r )

#define CNUM(a)    CNUM2(MKNAME(a,_re), MKNAME(a,_im))
#define CNUM2(a,b) (* (cnum_t*) & (num_t[2]) { a, b })

// --- num

int mad_num_sign (num_t x)
{
  return (x > 0) - (x < 0);  // -1, 0, 1
}

int mad_num_sign1 (num_t x)
{
  return ((int[]){ -1, 1 })[!signbit(x)]; // -1, 1: works for ±0, ±inf and ±NaN
}

num_t mad_num_fact (int n)
{
  int s = 1;

  if (n < 0) n = -n, s = n & 1 ? -s : s;

  return s*fact(n);
}

num_t mad_num_invfact (int n)
{
  int s = 1;

  if (n < 0) n = -n, s = n & 1 ? -s : s;

  return s/fact(n);
}

num_t mad_num_sinc (num_t x)
{
  return fabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x : sin(x)/x;
}

num_t mad_num_sinhc (num_t x)
{
  return fabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x : sinh(x)/x;
}

num_t mad_num_asinc (num_t x)
{
  return fabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x : asin(x)/x;
}

num_t mad_num_asinhc (num_t x)
{
  return fabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x : asinh(x)/x;
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

cnum_t mad_cnum_div (cnum_t x, cnum_t y)
{
// REFERENCES
//
// [1] Robert L. Smith. Algorithm 116: Complex division.  Commun. ACM,
//  5(8):435, 1962.
//
// [2] Michael Baudin and Robert L. Smith. "A robust complex division in
// Scilab," October 2012, available at http://arxiv.org/abs/1210.4539.

#define RBIG     (DBL_MAX/2)
#define RMIN     (DBL_MIN)
#define RMIN2    (0x1.0p-512)
#define RMINSCAL (0x1.0p+510)

  num_t x_re = creal(x), x_im = cimag(x);
  num_t y_re = creal(y), y_im = cimag(y);
  num_t r_re, r_im, denom, ratio;

  if (fabs(y_re) < fabs(y_im)) {
    /* avoid overflow/underflow issues when y_re and y_im are small */
    if (fabs(y_im) < RMIN2) {
      x_re = x_re * RMINSCAL;
      x_im = x_im * RMINSCAL;
      y_re = y_re * RMINSCAL;
      y_im = y_im * RMINSCAL;
    }
//    /* prevent overflow when arguments are near max representable */
//    else if ((fabs(y_im) >= RBIG) || (fabs(x_re) >= RBIG) || (fabs(x_im) >= RBIG)) {
//      x_re = x_re * 0.5;
//      x_im = x_im * 0.5;
//      y_re = y_re * 0.5;
//      y_im = y_im * 0.5;
//    }
    ratio = y_re / y_im;
    denom = (y_re * ratio) + y_im;
    if (fabs(ratio) > RMIN) {
      r_re = ((x_re * ratio) + x_im) / denom;
      r_im = ((x_im * ratio) - x_re) / denom;
    } else {
      r_re = ((y_re * (x_re / y_im)) + x_im) / denom;
      r_im = ((y_re * (x_im / y_im)) - x_re) / denom;
    }
  } else {
    /* avoid overflow issues when y_re and y_im are small */
    if (fabs(y_re) < RMIN2) {
      x_re = x_re * RMINSCAL;
      x_im = x_im * RMINSCAL;
      y_re = y_re * RMINSCAL;
      y_im = y_im * RMINSCAL;
    }
//   /* prevent overflow when arguments are near max representable */
//    else if ((fabs(y_re) >= RBIG) || (fabs(x_re) >= RBIG) || (fabs(x_im) >= RBIG)) {
//      x_re = x_re * 0.5;
//      x_im = x_im * 0.5;
//      y_re = y_re * 0.5;
//      y_im = y_im * 0.5;
//    }
    ratio = y_im / y_re;
    denom = (y_im * ratio) + y_re;
    if (fabs(ratio) > RMIN) {
      r_re = (x_re + (x_im * ratio)) / denom;
      r_im = (x_im - (x_re * ratio)) / denom;
    } else {
      r_re = (x_re + (y_im * (x_im / y_re))) / denom;
      r_im = (x_im - (y_im * (x_re / y_re))) / denom;
    }
  }
  return CNUM(r);

#undef RBIG
#undef RMIN
#undef RMIN2
#undef RMINSCAL
}

cnum_t mad_cnum_inv (cnum_t x)
{
  return mad_cnum_div(1, x);
}

cnum_t mad_cnum_sinc  (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x
                      : mad_cnum_div(csin(x), x);
}

cnum_t mad_cnum_sinhc (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x
                      : mad_cnum_div(csinh(x), x);
}

cnum_t mad_cnum_asinc  (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 + 0.1666666666666666666667*x*x
                      : mad_cnum_div(casin(x), x);
}

cnum_t mad_cnum_asinhc (cnum_t x)
{
  return cabs(x)<1e-4 ? 1 - 0.1666666666666666666667*x*x
                      : mad_cnum_div(casinh(x), x);
}

cnum_t mad_cnum_powi (cnum_t x, int n)
{
  cnum_t r = 1;
  if (n < 0) n = -n, x = mad_cnum_inv(x);
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

void mad_cnum_sinc_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_sinc  ( CNUM(x) ); }
void mad_cnum_sinhc_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_sinhc ( CNUM(x) ); }
void mad_cnum_asinc_r (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_asinc ( CNUM(x) ); }
void mad_cnum_asinhc_r(num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = mad_cnum_asinhc( CNUM(x) ); }

void mad_cnum_powi_r  (num_t x_re, num_t x_im, int n, cnum_t *r)
{ CHKR; *r = mad_cnum_powi( CNUM(x), n ); }

void mad_cnum_unit_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / cabs( CNUM(x) ); }

void mad_cnum_rect_r (num_t rho, num_t ang, cnum_t *r)
{ CHKR; *r = CNUM2( rho * cos(ang), rho * sin(ang) ); }

void mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM2( cabs(CNUM(x)), carg(CNUM(x)) ); }

void mad_cnum_div_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = mad_cnum_div(CNUM(x), CNUM(y)); }

void mad_cnum_inv_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = mad_cnum_div(1, CNUM(x)); }

void mad_cnum_mod_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; cnum_t cr = mad_cnum_div(CNUM(x), CNUM(y));
  *r = CNUM(x) - CNUM(y) * CNUM2(round(creal(cr)), round(cimag(cr))); }

void mad_cnum_pow_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = cpow( CNUM(x), CNUM(y) ); }

// --- Faddeeva function and variants from MIT --------------------------------o

#include "mad_erfw.h"

num_t  mad_num_wf     (num_t x) { return Faddeeva_w_im     (x); }
num_t  mad_num_erf    (num_t x) { return Faddeeva_erf_re   (x); }
num_t  mad_num_erfc   (num_t x) { return Faddeeva_erfc_re  (x); }
num_t  mad_num_erfi   (num_t x) { return Faddeeva_erfi_re  (x); }
num_t  mad_num_erfcx  (num_t x) { return Faddeeva_erfcx_re (x); }
num_t  mad_num_dawson (num_t x) { return Faddeeva_Dawson_re(x); }

cnum_t mad_cnum_wf    (cnum_t x, num_t relerr) { return Faddeeva_w     (x, relerr); }
cnum_t mad_cnum_erf   (cnum_t x, num_t relerr) { return Faddeeva_erf   (x, relerr); }
cnum_t mad_cnum_erfc  (cnum_t x, num_t relerr) { return Faddeeva_erfc  (x, relerr); }
cnum_t mad_cnum_erfi  (cnum_t x, num_t relerr) { return Faddeeva_erfi  (x, relerr); }
cnum_t mad_cnum_erfcx (cnum_t x, num_t relerr) { return Faddeeva_erfcx (x, relerr); }
cnum_t mad_cnum_dawson(cnum_t x, num_t relerr) { return Faddeeva_Dawson(x, relerr); }

void mad_cnum_wf_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_w (CNUM(x), relerr); }

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

// -- RNG XoShiRo256** --------------------------------------------------------o

#define N 4

static inline u64_t splitmix64 (u64_t x) {
  u64_t z = x         + 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline u64_t rotl(const u64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

u64_t mad_num_randi (prng_state_t *restrict st)
{
  u64_t *s = st->s;
  const u64_t r = rotl(s[1] * 5, 7) * 9;
  const u64_t t = s[1] << 17;

  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = rotl(s[3], 45);

  return r; // number within [0,ULLONG_MAX]
}

union numbit { u64_t u; num_t d; };

num_t mad_num_rand (prng_state_t *restrict st)
{
  u64_t x = mad_num_randi(st);
  const union numbit n = { .u = 0x3ffULL << 52 | x >> 12 };
  return n.d - 1; // [1.,2.) -> [0.,1.)
}

void mad_num_randseed (prng_state_t *restrict st, num_t seed)
{
  const union numbit n = { .d = seed };
  st->s[0] = splitmix64(n.u);
  for (int i=1; i < N; i++) st->s[i] = splitmix64(st->s[i-1]);
}

void mad_num_randjump (prng_state_t *restrict st)
{
  static const u64_t jump[N] = {
    0x180ec6d33cfd0abaULL, 0xd5a61266f0c9392cULL,
    0xa9582618e03fc9aaULL, 0x39abdc4529b1661cULL
  };

  u64_t *s = st->s;
  u64_t s0=0, s1=0, s2=0, s3=0;
  for(int i=0; i < N; i++)
    for(int b=0; b < 64; b++) {
      if (jump[i] & 1llu << b) {
        s0 ^= s[0];
        s1 ^= s[1];
        s2 ^= s[2];
        s3 ^= s[3];
      }
      mad_num_randi(st);
    }

  s[0] = s0;
  s[1] = s1;
  s[2] = s2;
  s[3] = s3;
}

#undef N

// -- RNG MADX ----------------------------------------------------------------o

#define MAX_RAND 1000000000
#define SCL_RAND 1e-9
#define  NR_RAND 55
#define  NJ_RAND 24
#define  ND_RAND 21

static void
irngen(xrng_state_t *rng)
{
  int j;
  for (int i = 0; i < NJ_RAND; i++) {
    if ((j = rng->s[i] - rng->s[i+NR_RAND-NJ_RAND]) < 0) j += MAX_RAND;
    rng->s[i] = j;
  }
  for (int i = NJ_RAND; i < NR_RAND; i++) {
    if ((j = rng->s[i] - rng->s[i-NJ_RAND]) < 0) j += MAX_RAND;
    rng->s[i] = j;
  }
  rng->n = 0;
}

void
mad_num_xrandseed (xrng_state_t *rng, u32_t seed)
{
  int k = 1, j = seed % MAX_RAND;
  rng->s[NR_RAND-1] = j;
  for (int i = 0; i < NR_RAND-1; i++) {
    int ii = (ND_RAND*(i+1)) % NR_RAND;
    rng->s[ii-1] = k;
    if ((k = j - k) < 0) k += MAX_RAND;
    j = rng->s[ii-1];
  }
  for (int i = 0; i < 3; i++) irngen(rng);
}

num_t
mad_num_xrand (xrng_state_t *rng) // [0.,1.)
{
  if (rng->n == NR_RAND) irngen(rng);
  return SCL_RAND * rng->s[rng->n++];
}

u32_t
mad_num_xrandi (xrng_state_t *rng) // [0,UINT_MAX]
{
  if (rng->n == NR_RAND) irngen(rng);
  return rng->s[rng->n++];
}

// -- OPENMP TEST -------------------------------------------------------------o

num_t mad_num_suminv(u64_t n)
{
  num_t s=0;

#ifdef _OPENMP
  #pragma omp parallel for reduction(+ : s)
#endif
  for (u64_t i=1; i <= n; ++i)
    s += 1.0/i;

  return s;
}
