/*
 o----------------------------------------------------------------------------o
 |
 | Number module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

#include <complex.h>
#include <assert.h>

#include "mad_num.h"

#define  num_t double
#define cnum_t double _Complex

// --- implementation --------------------------------------------------------o

#define CHKR  assert( r )

#define CNUM(a) (* (cnum_t*) & (num_t[2]) { (a##_re), (a##_im) })

// --- num

// --- cnum

num_t mad_cnum_abs  (num_t x_re, num_t x_im)  { return cabs( CNUM(x) ); }
num_t mad_cnum_arg  (num_t x_re, num_t x_im)  { return carg( CNUM(x) ); }

void mad_cnum_exp   (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = cexp   ( CNUM(x) ); }
void mad_cnum_log   (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = clog   ( CNUM(x) ); }
void mad_cnum_sqrt  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = csqrt  ( CNUM(x) ); }
void mad_cnum_proj  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = cproj  ( CNUM(x) ); }
void mad_cnum_sin   (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = csin   ( CNUM(x) ); }
void mad_cnum_cos   (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = ccos   ( CNUM(x) ); }
void mad_cnum_tan   (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = ctan   ( CNUM(x) ); }
void mad_cnum_sinh  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = csinh  ( CNUM(x) ); }
void mad_cnum_cosh  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = ccosh  ( CNUM(x) ); }
void mad_cnum_tanh  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = ctanh  ( CNUM(x) ); }
void mad_cnum_asin  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = casin  ( CNUM(x) ); }
void mad_cnum_acos  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = cacos  ( CNUM(x) ); }
void mad_cnum_atan  (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = catan  ( CNUM(x) ); }
void mad_cnum_asinh (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = casinh ( CNUM(x) ); }
void mad_cnum_acosh (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = cacosh ( CNUM(x) ); }
void mad_cnum_atanh (num_t x_re, num_t x_im, cnum_t *r)  { CHKR; *r = catanh ( CNUM(x) ); }

void mad_cnum_div (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = CNUM(x) / CNUM(y);  }

void mad_cnum_pow (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r)
{ CHKR; *r = cpow( CNUM(x), CNUM(y) ); }
