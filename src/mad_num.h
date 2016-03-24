#ifndef MAD_NUM_H
#define MAD_NUM_H

/*
 o----------------------------------------------------------------------------o
 |
 | Number module interface
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
  
  Purpose:
  - wrappers around functions of real and complex numbers for LuaJIT

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

num_t mad_cnum_abs_r   (num_t x_re, num_t x_im);
num_t mad_cnum_arg_r   (num_t x_re, num_t x_im);

void  mad_cnum_rect_r  (num_t  rho, num_t  ang, cnum_t *r);
void  mad_cnum_polar_r (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_exp_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_log_r   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_sqrt_r  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_proj_r  (num_t x_re, num_t x_im, cnum_t *r);

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

void  mad_cnum_ipow_r  (num_t x_re, num_t x_im, long long y, cnum_t *r);

void  mad_cnum_div_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_mod_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_pow_r   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);

// ---------------------------------------------------------------------------o

#endif // MAD_NUM_H