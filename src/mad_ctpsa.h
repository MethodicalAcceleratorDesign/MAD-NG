#ifndef MAD_CTPSA_H
#define MAD_CTPSA_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Complex Truncated Power Series Algebra module interface
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - provide a full feathered Generalized Complex TPSA package

  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad_mono.h"
#include "mad_desc.h"

// --- types -----------------------------------------------------------------o

struct tpsa;

typedef struct ctpsa ctpsa_t;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// --- interface -------------------------------------------------------------o

// ctors, dtor
ctpsa_t* mad_ctpsa_newd    (desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
ctpsa_t* mad_ctpsa_new     (const ctpsa_t *t, ord_t mo);
void     mad_ctpsa_del     (      ctpsa_t *t);

// introspection
desc_t*  mad_ctpsa_desc    (const ctpsa_t *t);
ord_t    mad_ctpsa_ord     (const ctpsa_t *t);
ord_t    mad_ctpsa_ordv    (const ctpsa_t *t1, const ctpsa_t *t2, ...);  // max order of all

// initialization
void     mad_ctpsa_copy    (const ctpsa_t *t, ctpsa_t *dst);
void     mad_ctpsa_clear   (      ctpsa_t *t);
void     mad_ctpsa_scalar  (      ctpsa_t *t, cnum_t v);
void     mad_ctpsa_scalar_r(      ctpsa_t *t, num_t v_re, num_t v_im); // without complex-by-value

// conversion
void     mad_ctpsa_real    (const ctpsa_t *t, struct tpsa *dst);
void     mad_ctpsa_imag    (const ctpsa_t *t, struct tpsa *dst);

// indexing / monomials
int      mad_ctpsa_mono    (const ctpsa_t *t, ssz_t n,       ord_t m_[n], idx_t i);
idx_t    mad_ctpsa_midx    (const ctpsa_t *t, ssz_t n, const ord_t m [n]);
idx_t    mad_ctpsa_midx_s  (const ctpsa_t *t, ssz_t n,       str_t s    ); // string mono "[0-9]*"
idx_t    mad_ctpsa_midx_sp (const ctpsa_t *t, ssz_t n, const int   m [n]); // sparse mono [(i,o)]

// accessors
cnum_t   mad_ctpsa_get0    (const ctpsa_t *t);
cnum_t   mad_ctpsa_geti    (const ctpsa_t *t, idx_t i);
cnum_t   mad_ctpsa_getm    (const ctpsa_t *t, ssz_t n, const ord_t m[n]);
cnum_t   mad_ctpsa_getm_s  (const ctpsa_t *t, ssz_t n,       str_t s   ); // string w orders in '0'-'9'
cnum_t   mad_ctpsa_getm_sp (const ctpsa_t *t, ssz_t n, const int   m[n]); // sparse mono [(i,o)]
void     mad_ctpsa_set0    (      ctpsa_t *t, /* i = 0 */                cnum_t a, cnum_t b); // a*x[0]+b
void     mad_ctpsa_seti    (      ctpsa_t *t, idx_t i,                   cnum_t a, cnum_t b); // a*x[i]+b
void     mad_ctpsa_setm    (      ctpsa_t *t, ssz_t n, const ord_t m[n], cnum_t a, cnum_t b); // a*x[m]+b
void     mad_ctpsa_setm_s  (      ctpsa_t *t, ssz_t n,       str_t s   , cnum_t a, cnum_t b); // a*x[m]+b
void     mad_ctpsa_setm_sp (      ctpsa_t *t, ssz_t n, const int   m[n], cnum_t a, cnum_t b); // a*x[m]+b

// accessors without complex-by-value
void     mad_ctpsa_get0_r   (const ctpsa_t *t, cnum_t *r);
void     mad_ctpsa_geti_r   (const ctpsa_t *t, idx_t i, cnum_t *r);
void     mad_ctpsa_getm_r   (const ctpsa_t *t, ssz_t n, const ord_t m[n], cnum_t *r);
void     mad_ctpsa_getm_s_r (const ctpsa_t *t, ssz_t n,       str_t s   , cnum_t *r);
void     mad_ctpsa_getm_sp_r(const ctpsa_t *t, ssz_t n, const int   m[n], cnum_t *r);
void     mad_ctpsa_set0_r   (      ctpsa_t *t, /* i = 0 */                num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_seti_r   (      ctpsa_t *t, idx_t i,                   num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setm_r   (      ctpsa_t *t, ssz_t n, const ord_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setm_s_r (      ctpsa_t *t, ssz_t n,       str_t s   , num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setm_sp_r(      ctpsa_t *t, ssz_t n, const int   m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im);

// operations
void     mad_ctpsa_abs     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_arg     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_conj    (const ctpsa_t *a, ctpsa_t *c);
cnum_t   mad_ctpsa_nrm1    (const ctpsa_t *a, const ctpsa_t *b_);
cnum_t   mad_ctpsa_nrm2    (const ctpsa_t *a, const ctpsa_t *b_);
void     mad_ctpsa_der     (const ctpsa_t *a, ctpsa_t *c, int var);  // TODO: check functions that rely on it
void     mad_ctpsa_mder    (const ctpsa_t *a, ctpsa_t *c, ssz_t n, const ord_t m[]);

void     mad_ctpsa_add     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_sub     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mul     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_div     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);

void     mad_ctpsa_acc     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c += v*a, aliasing OK
void     mad_ctpsa_scl     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v*a
void     mad_ctpsa_inv     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/a
void     mad_ctpsa_invsqrt (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/sqrt(a)

void     mad_ctpsa_sqrt    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_exp     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_log     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sin     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cos     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sinh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cosh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sincos  (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sincosh (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sinc    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sirx    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_corx    (const ctpsa_t *a, ctpsa_t *c);

void     mad_ctpsa_tan     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cot     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asin    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acos    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atan    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acot    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_tanh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_coth    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asinh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acosh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atanh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acoth   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_erf     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_ipow    (const ctpsa_t *a, ctpsa_t *c, int n);

// operations without complex-by-value
void     mad_ctpsa_nrm1_r   (const ctpsa_t *a, const ctpsa_t *b_, cnum_t *r);
void     mad_ctpsa_nrm2_r   (const ctpsa_t *a, const ctpsa_t *b_, cnum_t *r);
void     mad_ctpsa_acc_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_scl_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_inv_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_invsqrt_r(const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);

// high level functions
void     mad_ctpsa_axpb        (cnum_t a, const ctpsa_t *x,
                                cnum_t b, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axpbypc     (cnum_t a, const ctpsa_t *x,
                                cnum_t b, const ctpsa_t *y,
                                cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypb       (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbzpc    (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, const ctpsa_t *z,
                                cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbvwpc   (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, const ctpsa_t *v, const ctpsa_t *w,
                                cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_ax2pby2pcz2 (cnum_t a, const ctpsa_t *x,
                                cnum_t b, const ctpsa_t *y,
                                cnum_t c, const ctpsa_t *z, ctpsa_t *r); // aliasing OK

// high level functions without complex-by-value
void     mad_ctpsa_axpb_r       (num_t a_re, num_t a_im, const ctpsa_t *x,
                                 num_t b_re, num_t b_im, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axpbypc_r    (num_t a_re, num_t a_im, const ctpsa_t *x,
                                 num_t b_re, num_t b_im, const ctpsa_t *y,
                                 num_t c_re, num_t c_im, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypb_r      (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                 num_t b_re, num_t b_im, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbzpc_r   (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                 num_t b_re, num_t b_im, const ctpsa_t *z,
                                 num_t c_re, num_t c_im, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbvwpc_r  (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                 num_t b_re, num_t b_im, const ctpsa_t *v, const ctpsa_t *w,
                                 num_t c_re, num_t c_im, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_ax2pby2pcz2_r(num_t a_re, num_t a_im, const ctpsa_t *x,
                                 num_t b_re, num_t b_im, const ctpsa_t *y,
                                 num_t c_re, num_t c_im, const ctpsa_t *z, ctpsa_t *r); // aliasing OK

// to check for non-homogeneous maps & knobs
void     mad_ctpsa_poisson (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int n);  // TO CHECK n
void     mad_ctpsa_compose (int sa, const ctpsa_t *ma[], int sb, const ctpsa_t *mb[], int sc, ctpsa_t *mc[]);
void     mad_ctpsa_minv    (int sa, const ctpsa_t *ma[],                              int sc, ctpsa_t *mc[]);
void     mad_ctpsa_pminv   (int sa, const ctpsa_t *ma[],                              int sc, ctpsa_t *mc[], int row_select[]);

// I/O
void     mad_ctpsa_print    (const ctpsa_t *t, str_t name_, FILE *stream_);
ctpsa_t* mad_ctpsa_scan     (                               FILE *stream_); // TODO
desc_t*  mad_ctpsa_scan_hdr (                               FILE *stream_);
void     mad_ctpsa_scan_coef(      ctpsa_t *t,              FILE *stream_); // TODO
void     mad_ctpsa_debug    (const ctpsa_t *t);

#define  mad_ctpsa_ordv(...) mad_ctpsa_ordv(__VA_ARGS__,NULL)

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_H
