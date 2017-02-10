/*
 o-----------------------------------------------------------------------------o
 |
 | Number module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

// --- implementation ---------------------------------------------------------o

#define CHKR  assert( r )

#define CNUM2(a,b) (* (cnum_t*) & (num_t[2]) { a, b })
#define CNUM(a) CNUM2(MKNAME(a,_re), MKNAME(a,_im))

// --- cnum

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

void mad_cnum_unit_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / cabs( CNUM(x) ); }

void mad_cnum_rect_r (num_t rho, num_t ang, cnum_t *r)
{ CHKR; *r = CNUM2( rho * cos(ang), rho * sin(ang) ); }

void mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM2( cabs(CNUM(x)), carg(CNUM(x)) ); }

void mad_cnum_div_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / CNUM(y);  }

void mad_cnum_mod_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; cnum_t cr = CNUM(x) / CNUM(y);
  *r = CNUM(x) - CNUM(y) * CNUM2(round(creal(cr)), round(cimag(cr))); }

void mad_cnum_pow_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = cpow( CNUM(x), CNUM(y) ); }

void mad_cnum_ipow_r (num_t x_re, num_t x_im, long long y, cnum_t *r)
{ CHKR; cnum_t x = CNUM(x);
  *r = 1;
  if (y < 0) y = -y, x = 1/x;
  for (;;) {
    if (y & 1) *r *= x;
    if (y >>= 1) x *= x; else break;
  }
}

// --- Faddeeva function and variants from MIT

#include "Faddeeva.h"

void mad_cnum_erf_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erf (CNUM(x), relerr); }

void mad_cnum_erfc_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfc (CNUM(x), relerr); }

void mad_cnum_erfi_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfi (CNUM(x), relerr); }

void mad_cnum_erfw_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_w (CNUM(x), relerr); }

void mad_cnum_erfcx_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_erfcx (CNUM(x), relerr); }

void mad_cnum_dawson_r (num_t x_re, num_t x_im, num_t relerr, cnum_t *r)
{ CHKR; *r = Faddeeva_Dawson (CNUM(x), relerr); }
