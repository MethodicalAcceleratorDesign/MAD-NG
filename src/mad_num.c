/*
 o-----------------------------------------------------------------------------o
 |
 | Number module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015+
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

// --- num

num_t mad_num_asinh   (num_t x) { return asinh(x);  }
num_t mad_num_acosh   (num_t x) { return acosh(x);  }
num_t mad_num_atanh   (num_t x) { return atanh(x);  }

num_t mad_num_erf     (num_t x) { return erf(x);    }
num_t mad_num_tgamma  (num_t x) { return tgamma(x); }
num_t mad_num_lgamma  (num_t x) { return lgamma(x); }

// --- cnum

num_t mad_cnum_carg_r (num_t x_re, num_t x_im) { return carg( CNUM(x) ); }
num_t mad_cnum_norm_r (num_t x_re, num_t x_im) { return cabs( CNUM(x) ); }

void mad_cnum_proj_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cproj  ( CNUM(x) ); }

void mad_cnum_exp_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = cexp   ( CNUM(x) ); }
void mad_cnum_log_r   (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = clog   ( CNUM(x) ); }
void mad_cnum_sqrt_r  (num_t x_re, num_t x_im, cnum_t *r) { CHKR; *r = csqrt  ( CNUM(x) ); }

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

void mad_cnum_unit_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / cabs( CNUM(x) ); }

void mad_cnum_rect_r (num_t rho, num_t ang, cnum_t *r)
{ CHKR; *r = CNUM2( rho * cos(ang), rho * sin(ang) ); }

void mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r)
{ CHKR; *r = CNUM2( cabs(CNUM(x)), carg(CNUM(x)) ); }

void mad_cnum_div_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / CNUM(y);  }

void mad_cnum_mod_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; num_t _[1]; cnum_t cr = CNUM(x) / CNUM(y);
  *r = CNUM(y) * CNUM2(modf(creal(cr),_), modf(cimag(cr),_)); }

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
